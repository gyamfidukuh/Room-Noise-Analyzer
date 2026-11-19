%% ============================================================
%% ROOM ACOUSTICS ANALYSER
%% Measures RT60 per octave band (125 Hz - 8 kHz) + noise floor,
%% classifies the room, and writes a JSON result into a
%% OneDrive-synced folder for the Power Automate flow to pick up.
%% ============================================================

clear; clc; close all;

SHOW_DEBUG_PLOTS = false;   % set true only when manually debugging on a laptop

%% ============================================================
%% 1. ROOM IDENTIFICATION — must run first, before any audio I/O
%% ============================================================
scriptFolder = fileparts(mfilename('fullpath'));
configPath   = fullfile(scriptFolder, 'room_config.txt');
counterPath  = fullfile(scriptFolder, 'room_counter.txt');

if isfile(configPath)
    roomID = string(strtrim(fileread(configPath)));
else
    roomID = "";
end

if roomID ~= ""
    % --- EXISTING ROOM: reuse the ID already saved for this laptop ---
    fprintf("Loaded existing Room ID: %s\n", roomID);

else
    % --- NEW ROOM: auto-assign the next available number ---
    if isfile(counterPath)
        nextNum = str2double(strtrim(fileread(counterPath)));
        if isnan(nextNum) || nextNum < 1
            nextNum = 1;   % counter file corrupted/empty — restart safely
        end
    else
        nextNum = 1;       % first-ever run, no counter file yet
    end

    roomID = "Room" + string(sprintf('%02d', nextNum));

    fid = fopen(counterPath, 'w');
    if fid == -1
        error("Could not write room_counter.txt in %s", scriptFolder);
    end
    fprintf(fid, '%d', nextNum + 1);
    fclose(fid);

    fid2 = fopen(configPath, 'w');
    if fid2 == -1
        error("Could not write room_config.txt in %s", scriptFolder);
    end
    fprintf(fid2, '%s', roomID);
    fclose(fid2);

    fprintf("Registered NEW room: %s (saved to room_config.txt)\n", roomID);
end

resp = input(sprintf("Confirm measuring room '%s'? (Y/N): ", roomID), 's');
if ~strcmpi(resp, "Y")
    error("Measurement cancelled by user.");
end
fprintf("Measuring room: %s\n\n", roomID);

%% ============================================================
%% 2. CONFIGURATION
%% ============================================================
fs         = 48000;    % sample rate (Hz)
sweepDur   = 5;         % sweep duration (s)
fStart     = 20;        % sweep start freq (Hz)
fEnd       = 20000;     % sweep end freq (Hz)
tailDur    = 2;         % extra silence after sweep to capture reverb tail (s)
silenceDur = 5;         % noise floor recording duration (s)

octaveBands = [125, 250, 500, 1000, 2000, 4000, 8000]; % Hz, band centers

% NOTE: MATLAB struct field names cannot start with a digit, so these are
% prefixed with "RT60_" — this also matches your SharePoint column names
% (RT60_125Hz, RT60_250Hz, ...) so the JSON keys map straight across.
bandNames = ["RT60_125Hz","RT60_250Hz","RT60_500Hz","RT60_1kHz","RT60_2kHz","RT60_4kHz","RT60_8kHz"];

% Classification thresholds (per instructor spec)
RT60_GOOD_MAX   = 0.6;    % seconds — below this + quiet = good for speech
RT60_POOR_MIN   = 1.0;    % seconds — at/above this = too echoey regardless of noise
NOISE_FLOOR_MAX = -40;    % dB — at/above this = too noisy regardless of RT60 (tune to your mic)

% IMPORTANT: replace with your actual OneDrive-synced path
oneDriveFolderPath = "C:\Users\USER\OneDrive - Kwame Nkrumah Uni. of Sci. and Tech\AcousticsData";


if ~isfolder(oneDriveFolderPath)
    error("OneDrive output folder not found: %s\nCheck the path and that OneDrive is synced.", oneDriveFolderPath);
end

%% ============================================================
%% 3. GENERATE TEST SIGNAL — logarithmic sweep (20 Hz to 20 kHz)
%% ============================================================
fprintf("Generating logarithmic sweep...\n");

t = (0:1/fs:sweepDur - 1/fs)';
sweepSignal = chirp(t, fStart, sweepDur, fEnd, 'logarithmic');

% Pad with silence so the room's reverb tail is fully captured
tailSamples     = tailDur * fs;
playbackSignal  = [sweepSignal; zeros(tailSamples, 1)];

%% ============================================================
%% 4. PLAY + RECORD SIMULTANEOUSLY
%% ============================================================
fprintf("Playing sweep and recording room response...\n");

recObj = audiorecorder(fs, 16, 1);
disp('Measurement starting... please be quiet!');

record(recObj);
pause(0.2);                          % start mic slightly early
sound(playbackSignal, fs);
pause(length(playbackSignal)/fs + 0.5);   % wait for full playback + tail
stop(recObj);

recordedSignal = getaudiodata(recObj);
disp('Sweep measurement complete.');

if SHOW_DEBUG_PLOTS
    figure; plot(recordedSignal); title('Raw Room Response');
end

%% ============================================================
%% 5. RECOVER IMPULSE RESPONSE (FFT-based deconvolution)
%% ============================================================
fprintf("Deconvolving to recover impulse response...\n");

paddedSweep = [sweepSignal; zeros(length(recordedSignal) - length(sweepSignal), 1)];

NFFT = 2^nextpow2(length(recordedSignal));
fft_recorded = fft(recordedSignal, NFFT);
fft_sweep    = fft(paddedSweep, NFFT);

fft_impulse = fft_recorded ./ (fft_sweep + eps);
impulse_response_full = real(ifft(fft_impulse, NFFT));

decaySamples = tailDur * fs;                    % keep the reverb-tail window
impulse_response = impulse_response_full(1:decaySamples);
impulse_response = impulse_response / max(abs(impulse_response));   % normalize to 0 dBFS

if SHOW_DEBUG_PLOTS
    figure;
    plot((0:length(impulse_response)-1)/fs, impulse_response);
    title('Deconvolved Room Impulse Response');
    xlabel('Time (seconds)'); ylabel('Amplitude');
end

% Trim to the main peak onward (discard pre-ringing/leading silence)
[~, peakIdx] = max(abs(impulse_response));
truncated_impulse = impulse_response(peakIdx:end);

%% ============================================================
%% 6. RT60 PER OCTAVE BAND (125 Hz - 8 kHz, Schroeder integration)
%% ============================================================
fprintf("Computing RT60 per octave band...\n");

rt60ByBand = struct();

for i = 1:length(octaveBands)
    fc = octaveBands(i);

    % 1-octave-wide bandpass filter centered at fc
    octFilt = octaveFilter(fc, 'SampleRate', fs, 'FilterOrder', 6);
    bandSignal = octFilt(truncated_impulse);

    % Schroeder backward integration
    schroeder_curve = flipud(cumsum(flipud(bandSignal.^2)));
    decay_dB = 10 * log10(schroeder_curve / max(schroeder_curve) + eps);
    t_decay  = (0:length(decay_dB)-1)' / fs;

    % T20 method: slope between -5 dB and -25 dB, extrapolated to -60 dB
    idx_minus5  = find(decay_dB <= -5,  1, 'first');
    idx_minus25 = find(decay_dB <= -25, 1, 'first');

    if isempty(idx_minus5) || isempty(idx_minus25) || idx_minus25 <= idx_minus5
        rt60Val = NaN;   % not enough clean decay range captured in this band
    else
        rt60Val = (t_decay(idx_minus25) - t_decay(idx_minus5)) * 3;  % T20 -> RT60
    end

    rt60ByBand.(bandNames(i)) = round(rt60Val, 3);
    fprintf("  %-12s : RT60 = %.3f s\n", bandNames(i), rt60Val);

    if SHOW_DEBUG_PLOTS && fc == 1000
        figure; plot(t_decay, decay_dB, 'LineWidth', 1.5); grid on;
        title('Logarithmic Decay Curve (1 kHz Octave Band)');
        xlabel('Time (seconds)'); ylabel('Energy (dB)'); ylim([-60 0]);
    end
end

rt60Values  = struct2array(rt60ByBand);
rt60Average = round(mean(rt60Values(~isnan(rt60Values))), 3);
fprintf("  Average RT60 = %.3f s\n\n", rt60Average);

%% ============================================================
%% 7. NOISE FLOOR — record silence, take RMS level in dB
%% ============================================================
fprintf("Recording %d s of silence for noise floor...\n", silenceDur);
fprintf("(Please keep the room quiet)\n");

noiseRecorder = audiorecorder(fs, 16, 1);
record(noiseRecorder);
pause(silenceDur);
stop(noiseRecorder);

noiseData    = getaudiodata(noiseRecorder);
noiseRMS     = sqrt(mean(noiseData.^2));
noiseFloorDB = round(20 * log10(noiseRMS + eps), 1);

fprintf("Noise floor = %.1f dB\n\n", noiseFloorDB);

%% ============================================================
%% 8. CLASSIFY THE ROOM (RT60 + noise floor together)
%% ============================================================
if rt60Average < RT60_GOOD_MAX && noiseFloorDB < NOISE_FLOOR_MAX
    verdict = "Good";
elseif rt60Average >= RT60_POOR_MIN || noiseFloorDB >= NOISE_FLOOR_MAX
    verdict = "Poor";
else
    verdict = "Marginal";
end

fprintf("Verdict: %s\n\n", verdict);

%% ============================================================
%% 8.5 GET GPS COORDINATES (Windows Location Services via .NET)
%% ============================================================
fprintf("Getting GPS coordinates from Windows Location Services...\n");

gpsTimeoutSec = 10;   % how long to wait for a location fix before giving up
latitude  = [];
longitude = [];

try
    NET.addAssembly('System.Device');
    watcher = System.Device.Location.GeoCoordinateWatcher();
    watcher.Start();

    tStart = tic;
    while watcher.Status ~= System.Device.Location.GeoPositionStatus.Ready ...
            && toc(tStart) < gpsTimeoutSec
        pause(0.5);
    end

    if watcher.Status == System.Device.Location.GeoPositionStatus.Ready
        coord = watcher.Position.Location;
        latitude  = double(coord.Latitude);
        longitude = double(coord.Longitude);
        fprintf("GPS fix acquired: lat=%.6f, long=%.6f\n\n", latitude, longitude);
    else
        fprintf("No GPS fix within %d s — leaving coordinates empty.\n\n", gpsTimeoutSec);
    end

    watcher.Stop();
catch ME
    fprintf("GPS lookup failed (%s) — leaving coordinates empty.\n\n", ME.message);
    latitude  = [];
    longitude = [];
end

%% ============================================================
%% 9. WRITE JSON OUTPUT (field names match the SharePoint/Flow schema)
%% ============================================================
fprintf("Writing result JSON...\n");

result = struct();
result.roomId       = roomID;
result.timestamp    = string(datetime('now', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
result.rt60         = rt60ByBand;      % per-band struct: 125Hz, 250Hz, ... 8kHz
result.rt60Average  = rt60Average;
result.noiseFloor    = noiseFloorDB;
result.verdict       = verdict;

% GPS — populated above if Windows Location Services returned a fix,
% otherwise stays empty (Parse JSON/SharePoint already handle this as optional)
result.latitude  = latitude;
result.longitude = longitude;

jsonText = jsonencode(result, 'PrettyPrint', true);

fileTimestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
outFilename   = fullfile(oneDriveFolderPath, roomID + "_" + fileTimestamp + ".json");

fid = fopen(outFilename, 'w');
if fid == -1
    error("Could not open output file for writing: %s", outFilename);
end
fwrite(fid, jsonText, 'char');
fclose(fid);

fprintf("Done. Result written to:\n  %s\n", outFilename);