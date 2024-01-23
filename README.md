# Ultrasound C-Scan Processing Application

This application is designed for processing and visualizing ultrasound C-scan data. It offers functionalities to load, save, trim, and add noise to C-scan data, as well as visualize it in various forms including orthoslices.

## Features

- **Data Loading**: Load ultrasound C-scan data from various formats including `.txt`, `.csv`, `.npy`, and `.bin`.
- **Data Saving**: Save processed C-scan data back to the disk.
- **Visualization**: View orthoslices of the C-scan data and manipulate visual aspects like surface alignment and noise addition.
- **Data Trimming**: Trim the dataset to focus on areas of interest.
- **Surface Calculation**: Determine the front surface of the C-scan data for further analysis.
- **Noise Addition**: Add Gaussian noise to the C-scan data for analysis under different conditions.

## Requirements

- Qt 5.x
- C++ compiler compatible with C++11 or later

## Building and Running the Application

1. Clone the repository:
git clone https://github.com/your-repository/ultrasound-cscan-processing.git

markdown
Copy code
2. Open the project in Qt Creator or your preferred Qt development environment.
3. Build the project using the build tools provided by the IDE.
4. Run the application from the IDE or the generated executable.

## Usage

1. **Load Data**: Use the 'Load data' button to load C-scan data from supported file formats.
2. **View Data**: Use the 'Orthoslice' button to visualize the data in orthogonal slices.
3. **Trim Data**: Use the 'Trim the dataset' button to focus on specific areas within the C-scan data.
4. **Add Noise**: Use the 'Add Noise' button to add Gaussian noise to the dataset.
5. **Save Data**: Use the 'Save' button to save the processed data.
