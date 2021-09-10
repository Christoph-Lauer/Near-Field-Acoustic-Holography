% @title        FFT based near-field acoustic holography algorithm
% @file         NearFielAcousticHoloraphy.m
% @short        An 2d fft based implementation of an near fiel acoustic holography alogrihm
% @version      0.1
% @date         10. September 2021
% @copyright    All rights reserved by Christoph Lauer
% @author       Christoph Lauer
% @contributors persons
% @client       company
% @language     MATLAB R2021a (Octave) 
% @packages     none
% @param        none
% @return       none
% @notes        you dont need the original CSV file from the LOUD array, you can use the converted wave file
% Toolbox       none
% @todo         movie genrator
% @copyright    Christoph Lauer
% @license      GPL 3.0
% @brief        This matlab file implements an NAH algorithm and generates an image sequence
% @contact      christophlauer@me.com
% @webpage      https://christoph-lauer.github.io
% 
% We implement here an acoustic near field algorithm which shift the z-axis in the frequency domain. The input data came from the LOUD
% microphone array where the input data is public available, at least in the waback internet archive. The raw data for the array is 
% available as CSV files which are converted in multi channel wave files which is re-opened. Once the conversion is completed we can 
% only use the wave file. The geometry data for the array is also available and is used in the getArrayFromRawInput function to extract
% one spund-preassure-level (SPL) array which is aligned in the correct order of the microphone channels corresponding to the array coordinates.
% This function returns one X-Y array with the SPL at for ine sample time. The recorded singal is 48000 samples long with an samplerate of 16000 1/s
% and a length of 3 seconds. This means that the getArrayFromRawInput delivers a two dimensional marirx from the array for each time sample point.
% The main algorithm lops over the time samples and extracts a 3D matrix for each sample point which is used to generate the slice plot which is
% stored on disc. So theoretical it is possible to generate a time lapse image sequence from 48000 images. The main algorith transforms the 2D array in the
% frequency domain, and shifts the signal in the z-achsis for each sample point. The LOUD array has a dimension of 60 x 17 = 1020 microphones.
% The 2D gaussian window function can be used as spatial low pass in the frequency domain. The NAH-algorithm generates for each time sample from the array
% an 3D hologram from the shifted z-axis.


%% 0.0) CONSTANTS
sampleRate = 16000;
global numMicrophones;
global coordinates;
global arrayRawData;
NAH = [];

%% 1.0) CONVERT THE CSV TO A WAVE FILE
% NOTE --> this has to be done only once a time
% 1.1) parse CSV File
%csv = readtable('data/R7-chirp-1024.csv');
%arrayRawData = table2array(csv);
%sound(arrayRawData(:,1), sampleRate)
% 1.2) save wave file
%audiowrite('data/samples.wav', arrayRawData, sampleRate);


%% 2.0) OPEN THE WAVE FILE
% 2.1) prepare the raw audio data
[arrayRawData sampleRate] = audioread('data/samples.wav');
sound(arrayRawData(:,510), sampleRate);
numSamples = length(arrayRawData);
arrayRawData(:,[509,510,511,512]) = []; % remove the special channels 509...512
numMicrophones = width(arrayRawData);
% 2.2) open the microphone positions file
coordinates = dlmread('data/coordinates.text');


%% 3.0) NAH, FFT-BASED ALGORITHM
FILTER = gaussianWindow2D(60,17);
for sample = 1:600 % loop over the samples of the array (48000 = 3s possible)
    SPL = getArrayFromRawInput(sample);
    count = 0;
    for shiftZ = 0:0.005:0.5 % loop over the z shife                   % SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT 
        frequencyDomain = fft2(SPL);                                   % HIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT
        %frequencyDomain = frequencyDomain .* FILTER;                  % IFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT S
        frequencyDomain = frequencyDomain .* exp(-i*2*pi*shiftZ);      % FT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SH
        timeDomain = ifft2(frequencyDomain);                           % T SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHI
        count = count + 1;                                             %  SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIF
        NAH(count,:,:) = real(timeDomain);                             % SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT SHIFT
    end
    % slice plot
    slice(NAH, 10:20:50, 1:length(NAH)/3:length(NAH), 8)
    colormap(hsv)
    title('Near-Field Acoustic Holography (NAH)')
    xlabel('Microphone Array X') 
    ylabel('Distance Z') 
    zlabel('Microphone Array Y') 
    % save image to disk
    ax = gca;
    filename = [sprintf('slice%03d',sample) '.png'];
    exportgraphics(ax, filename)
end

%% SAMPLES FROM GEOMETRY 
% this function tries to get back an input array from the microphones array corresponding the microphone array COORDIANTES
function [imageArray] = getArrayFromRawInput(sampleNumber)
    global numMicrophones;
    global coordinates;
    global arrayRawData;
    imageArray = zeros(60,17);            % pre allocate the image array
    for m = 1:numMicrophones              % loop over the microphone channels
        posX = coordinates(m,2);          % read the X coordinate
        posY = coordinates(m,3);          % read the Y coordinate
        indexX = round(posX / 0.03) + 28; % get the X index
        indexY = round(posY / 0.03) + 1;  % get the Y index
        imageArray(indexX, indexY) = arrayRawData(sampleNumber, m);
    end
end

%% (Note: not really used) 2D GAUSSIAN WINDOW
function w = gaussianWindow2D(N,M)
    wc=window('gausswin',N);
    wr=window('gausswin',M);
    [maskr,maskc]=meshgrid(wr,wc);
    w=maskr.*maskc;
end