# Audio Transmission System using MATLAB

This repository contains MATLAB scripts for implementing an audio transmission system utilizing Amplitude Modulation (AM) and Binary Phase Shift Keying (BPSK). The provided scripts demonstrate the complete process of reading audio signals, filtering, quantization, modulation, transmission, and demodulation.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Prerequisites](#prerequisites)
3. [Project Structure](#project-structure)
4. [Usage](#usage)
5. [Scripts Explanation](#scripts-explanation)

## Getting Started

Clone the repository to your local machine using:

```bash
git clone https://github.com/yourusername/audio-transmission-system.git
```

Navigate to the project directory:

```bash
cd audio-transmission-system
```

## Prerequisites

Ensure you have MATLAB installed on your machine. This project was developed and tested with MATLAB R2020a, but any recent version should work.

You will also need audio files named `test.wav` and `test1.wav` in the project directory. These files are used as input signals for the modulation and transmission processes.

## Project Structure

The repository contains the following scripts:

- `amplitude_modulation.m`: Script for amplitude modulation, including filtering and quantization.
- `bpsk_modulation.m`: Script for BPSK modulation and demodulation.

## Scripts Explanation

### `amplitude_modulation.m`

This script performs the following tasks:

1. **Reading Audio Signal**: Reads the audio file `test1.wav`.
2. **Filtering**: Applies a Butterworth low-pass filter to remove high-frequency noise.
3. **Quantization**: Quantizes the filtered signal to reduce the number of bits required for transmission.
4. **Pulse Code Modulation (PCM)**: Converts the quantized signal to a binary format.
5. **BPSK Modulation**: Modulates the binary signal using Binary Phase Shift Keying.
6. **Transmission and Reception**: Adds noise to the transmitted signal and demodulates it to recover the original audio.

### `bpsk_modulation.m`

This script performs the following tasks:

1. **Reading Audio Signals**: Reads two audio files `test.wav` and `test1.wav`.
2. **Filtering**: Applies a Butterworth low-pass filter to both signals.
3. **Carrier Signals**: Generates carrier signals for modulation.
4. **Single Sideband (SSB) Modulation**: Modulates the filtered signals using Single Sideband modulation.
5. **Multiplexing**: Combines the modulated signals for transmission.
6. **Demodulation**: Demodulates the received signal to recover the original audio signals.

## Example Output

The scripts plot various signals at different stages of the modulation and demodulation process, including:

- Original Message Signal
- Quantized Signal
- PCM Plot
- BPSK Signal
- Demodulated Signal

These plots help visualize the transformations that the signal undergoes during the transmission process.
