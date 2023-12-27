%%
% Ionoacoustics simulation (LMU)
% by MM
% last update: 30 Sept 2023

% ---------
%addpath('/Users/mariamaxouti/desktop/k-Wave_Code/k-wave-toolbox-version-1/k-Wave');
%addpath('/Users/mariamaxouti/desktop/k-Wave_Code/k-wave-toolbox-version-1/k-Wave/kWaveArray_alpha_0');

addpath('/rds/general/user/mm2321/home/k-wave-toolbox-version-1.4/k-Wave');

%addpath('\Users\herac\OneDrive - University College London\Documents\Downloads\k-wave-toolbox-version-1.3\k-Wave');
%addpath('\Users\herac\OneDrive - University College London\Documents\Downloads\k-wave-toolbox-version-1.3\k-Wave\kWaveArray_alpha_0');

%addpath('/home/mm2321/Desktop/change_size_on_axis/k-wave-toolbox-version-1/k-Wave');
%addpath('/home/mm2321/Desktop/change_size_on_axis/k-wave-toolbox-version-1/k-Wave/kWaveArray_alpha_0');

%addpath('/Users/mariamaxouti/Desktop/k-Wave_Code/k-wave-toolbox-version-1/vol3d');

tic

% clears the variables from memory
clearvars; 

%% ========================================================================
% DEFINE LITERALS (all the simulation-specific numbers)
% =========================================================================

% Computational parameters
GRID_SPACING = 1e-4;        % [m] i.e. 0.1 mm

COMPUTATIONAL_GRID_SIZE = 220;   
PML_SIZE = 20;
CFL = 0.3;

% Reconstruction type:      Iterative Time-Reversal
NUMBER_OF_TR_ITERATIONS = 3;

% Material properties
GRUENEISEN = 0.11;

% Acoustic sensor
single_element_transducer = false;
SINGLE_ELEMENT_DIAMETER = 1.3e-2;   % [m], i.e. 13 mm, V303 SU

array_transducer = false;           % shape: hemispherical 
%ARRAY_ELEMENT_DIAMETER = 3e-4;     % [m]
%NUMBER_OF_ELEMENTS = 400;  
ARRAY_ELEMENT_DIAMETER = 1e-3;      % [m]
NUMBER_OF_ELEMENTS = 186;  
ARRAY_RADIUS = 8e-3;                % [m]
ARRAY_Z_OFFSET = 0.0015;            % [m]

matrix_array = true;                % Verasonics matrix transducer
column_array = false;               % RC6gV Row-Column Array transducer
linear_array = false;               % 9L-D transducer

on_Axis = true;
Side = false;

% Define the frequency response of the sensor elements
if matrix_array
    center_freq = 3.5e6;      % [Hz]
    %bandwidth = 60;          % [%]
    bandwidth_range = 0.5e6;  % [MHz]
    bandwidth = ((2 * bandwidth_range) / center_freq) * 100;
end
if column_array
    center_freq = 6e6;    % [Hz]
    bandwidth = 100;      % [%]
end
if linear_array
    center_freq = 5.3e6;  % [Hz]
    bandwidth = 75;       % [%]
end

% options: true, false
Kapton_window = true;     
SciFiStations = false;    
SensorResponse = true;    
Electronic_noise = true;

if Kapton_window
    if SciFiStations
        add = 'with_window_with_SciFi';
    else
        add = 'with_window_without_SciFi';
    end
else 
    if SciFiStations
        add = 'without_window_with_SciFi';
    else
        add = 'without_window_without_SciFi';
    end
end

if SensorResponse
        add_Response = '_ResponseOn';
    else
        add_Response ='';
end
if Electronic_noise
        %add_ElectronicNoise = '_ElectronicNoise_';
        add_ElectronicNoise = ['_ElectronicNoise_' num2str(center_freq/1e6) 'MeV'];
    else
        add_ElectronicNoise ='';
end

%% ========================================================================
% SET UP THE SIMULATION AND SOURCE
% =========================================================================

% defines number of grid points in the simulation, taking into account the PML absorbing layer
%Nx = COMPUTATIONAL_GRID_SIZE-PML_SIZE;
%Ny = COMPUTATIONAL_GRID_SIZE-PML_SIZE;
%Nz = COMPUTATIONAL_GRID_SIZE-PML_SIZE;

% loads energy distribution from G4 simulation [J]
load energy_data2_cpp.mat
e0_protons = energy_data2_cpp;  % [J]

load electrons.mat
e0_electrons = energy_data2_cpp;  % [J]

%e0 = smooth(e0, true);
e0 = e0_protons + e0_electrons;

% define the grid spacing
dx = GRID_SPACING;     % grid point spacing in the x direction [m]
dy = GRID_SPACING;     % grid point spacing in the y direction [m]
dz = GRID_SPACING;     % grid point spacing in the z direction [m]

% Grid size 
Nx = 155;
Ny = 155;
Nz = 155;

% Adjust the input data (assuming you have a variable named 'inputData')
e0 = e0(50:150, 50:150, 50:150);

% Desired final size
final_size = [Nx, Ny, Nz];
original_data = e0; 

% Calculate the padding size for each dimension
padding_x = final_size(1) - size(original_data, 1);
padding_y = final_size(2) - size(original_data, 2);
padding_z = final_size(3) - size(original_data, 3);

% Make sure the padding sizes are non-negative
padding_x = max(padding_x, 0);
padding_y = max(padding_y, 0);
padding_z = max(padding_z, 0);

% Create an empty array of the final size
padded_data = zeros(final_size);

% Calculate the starting index for copying the original data to be centered on the left side
start_x = (final_size(1) - size(original_data, 1)) / 2 + 1;
start_y = (final_size(2) - size(original_data, 2)) / 2 + 1;
start_z = 1;

% Copy the original data into the left-centered region of the padded array
padded_data(start_x:start_x+size(original_data, 1)-1, start_y:start_y+size(original_data, 2)-1, start_z:start_z+size(original_data, 3)-1) = original_data;

% Define the source energy as the new padded data
e0 = padded_data;

% create the grid object
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% divide energy by voxel volume to get energy density
voxel_volume = dx*dy*dz; % m^3
e0_density = e0/voxel_volume;         
% multiply by the dimensionless Grueneisen parameter to get initial
% pressure distribution and assign to the source structure 
source.p0 = GRUENEISEN * e0_density;

% if necessary smooth the data,  eg. to reduce the effect of using too few
% particles in Geant4
%source.p0 = smooth(source.p0, true);

p0 = source.p0;

% Swap X and Z axes
%p0 = permute(p0, [3, 2, 1]);

% Plot the 3D pressure profile
%fig = gcf;

% Enable interactivity - Data cursor
%datacursormode on;
%vol3d('cdata', p0, 'texture', '3D');

% Set labels and title
%xlabel('X [m]');
%ylabel('Y [m]');
%zlabel('Z [m]');
%title('Pressure Profile');

% Adjust view and aspect ratio
%view(3);
%view(-90,0);
%daspect([kgrid.dx, kgrid.dy, kgrid.dz]);  % Adjust aspect ratio based on grid spacing
%rotate3d on;

% Save the vol3d plot as a PNG file
%saveas(gcf, 'vol3d_plot.png', 'png');

%%%% integrating radially the initial energy distribution %%%%
% empty array to store values
int_x_e = zeros(1,Nx);
int_y_e = zeros(1,Ny);
int_z_e = zeros(1,Nz);

plane_xy_int = zeros(Nz);
plane_xz_int = zeros(Ny);
plane_yz_int = zeros(Nx);

% integrate along each plane
% and store total energy in each plane 
for e = 1:Nx
    
    e0_MeV = e0*6.242e12;   % [MeV]

    % pressure in each plane = 1 voxel
    plane_xy = squeeze(e0_MeV(:,:,e));
    plane_xz = squeeze(e0_MeV(:,e,:));
    plane_yz = squeeze(e0_MeV(e,:,:));
    
    % adding the pressure in all voxels
    combined_xy = sum(plane_xy, 'all'); 
    combined_xz = sum(plane_xz, 'all'); 
    combined_yz = sum(plane_yz, 'all'); 
    % store into arrays
    int_x_e(e) = combined_yz;
    int_y_e(e) = combined_xz;
    int_z_e(e) = combined_xy;

    % adds the images
    plane_xy_int = plane_xy_int + plane_xy;
    plane_xz_int = plane_xz_int + plane_xz;
    plane_yz_int = plane_yz_int + plane_yz;
end

filename = sprintf('int_e0_z_source.csv');   % [MeV]
%writematrix(int_z_e, filename); 

data = int_x_e;
halfMax = (min(data) + max(data))/2;  % half max value
% find where the data first drops below half the max
index1 = find(data >= halfMax, 1, 'first');
% find where the data last rises above half the max
index2 = find(data >= halfMax, 1, 'last');
fwhm_x = index2-index1 + 1;    % FWHM in indexes [voxels]
fwhm_x = fwhm_x*GRID_SPACING;  % [m]
fwhm_x = fwhm_x*1000;          % [mm]

data = int_y_e;
halfMax = (min(data) + max(data))/2; % half max value
% find where the data first drops below half the max
index1 = find(data >= halfMax, 1, 'first');
% find where the data last rises above half the max
index2 = find(data >= halfMax, 1, 'last');
fwhm_y = index2-index1 + 1;    % FWHM in indexes [voxels]
fwhm_y = fwhm_y*GRID_SPACING;  % [m]
fwhm_y = fwhm_y*1000;          % [mm]

x_vec = (0:Nx-1)*dx; % set x_vec to start from zero
y_vec = (0:Ny-1)*dy; % set y_vec to start from zero

% find the maximum energy
maxValue = max([max(plane_xy_int), max(plane_xz_int), max(plane_yz_int)]);

figure;
%sgtitle('Integrated Initial Energy Distribution')
subplot(2,2,1);
%imagesc(plane_xy_int);
imagesc(y_vec*1e3,x_vec*1e3,plane_xy_int); % plot in [mm]
title('x-y plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
caxis([0, maxValue]);
%set(gca,'ColorScale','log');
c.Label.String = 'Energy [MeV]';

subplot(2,2,2);
%imagesc(plane_xz_int);
imagesc(y_vec*1e3,x_vec*1e3,plane_xz_int); % plot in [mm]
title('x-z plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
caxis([0, maxValue]);
%set(gca,'ColorScale','log');
c.Label.String = 'Energy [MeV]';

subplot(2,2,3);
%imagesc(plane_yz_int);
imagesc(y_vec*1e3,x_vec*1e3,plane_yz_int); % plot in [mm]
title('y-z plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
caxis([0, maxValue]);
%set(gca,'ColorScale','log');
c.Label.String = 'Energy [MeV]';
f = gcf;
exportgraphics(f, ['integrated_energy_distribution_source', '.png'], 'Resolution', 300);


%%%% integrating radially the initial pressure distribution %%%%
% empty array to store values
int_x_p = zeros(1,Nx);
int_y_p = zeros(1,Ny);
int_z_p = zeros(1,Nz);

plane_xy_int = zeros(Nz);
plane_xz_int = zeros(Ny);
plane_yz_int = zeros(Nx);

% integrate along each plane
% and store total pressure in each plane - along z-axis
for p = 1:Nx
    
    % pressure in each plane = 1 voxel
    plane_xy = squeeze(p0(:,:,p));
    plane_xz = squeeze(p0(:,p,:));
    plane_yz = squeeze(p0(p,:,:));
    
    % adding the pressure in all voxels
    combined_xy = sum(plane_xy, 'all'); 
    combined_xz = sum(plane_xz, 'all'); 
    combined_yz = sum(plane_yz, 'all'); 
    % store into arrays
    int_x_p(p) = combined_yz;
    int_y_p(p) = combined_xz;
    int_z_p(p) = combined_xy;

    % adds the images
    plane_xy_int = plane_xy_int + plane_xy;
    plane_xz_int = plane_xz_int + plane_xz;
    plane_yz_int = plane_yz_int + plane_yz;
end

filename = sprintf('int_p0_x_source.csv');
writematrix(int_x_p, filename); 
filename = sprintf('int_p0_y_source.csv');
writematrix(int_y_p, filename); 
filename = sprintf('int_p0_z_source.csv');
writematrix(int_z_p, filename); 

% find the maximum pressure
maxValue = max([max(plane_xy_int), max(plane_xz_int), max(plane_yz_int)]);

figure;
%sgtitle('Integrated Initial Pressure Distribution')
subplot(2,2,1);
%imagesc(plane_xy_int);
imagesc(y_vec*1e3,x_vec*1e3,plane_xy_int); % plot in [mm]
title('x-y plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
caxis([0, maxValue]);
c.Label.String = 'Pressure [Pa]';

subplot(2,2,2);
%imagesc(plane_xz_int);
imagesc(y_vec*1e3,x_vec*1e3,plane_xz_int); % plot in [mm]
title('x-z plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
caxis([0, maxValue]);
c.Label.String = 'Pressure [Pa]';

subplot(2,2,3);
%imagesc(plane_yz_int);
imagesc(y_vec*1e3,x_vec*1e3,plane_yz_int); % plot in [mm]
title('y-z plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
caxis([0, maxValue]);
c.Label.String = 'Pressure [Pa]';
f = gcf;
exportgraphics(f, ['integrated_pressure_distribution_source', '.png'], 'Resolution', 300);

% for difference image later
plane_xy_int_source = plane_xy_int;
plane_xz_int_source = plane_xz_int;
plane_yz_int_source = plane_yz_int;

%% ========================================================================
% DEFINE KAPTON FOIL
% =========================================================================

WATER_DENSITY = 997;        % [kg/m^3]
SOUND_SPEED_WATER = 1500;   % [m/s]
medium.sound_speed = SOUND_SPEED_WATER*ones(Nx, Ny, Nz);
medium.density = WATER_DENSITY*ones(Nx, Ny, Nz);

if Kapton_window

    % define the kapton foil
    kapton_x_pos = Nx/2;
    kapton_y_pos = Ny/2;
    kapton_z_offset = 65;   % [grid points]
    kapton_z_pos = Nz/2 - kapton_z_offset;
    Kapton_foil = zeros(Nx, Ny, Nz);
    
    kapton_thickness = 1;          % [grid points], 1 grid = 0.1 mm
    kapton_x = 150;                % [grid points], 15 mm
    kapton_y = 150;                % [grid points], 15 mm
    kapton_z = kapton_thickness;   % [grid points], 0.1 mm
     
    Kapton_foil(kapton_x_pos - kapton_x/2 : kapton_x_pos + kapton_x/2, ...
        kapton_y_pos - kapton_y/2 : kapton_y_pos+ kapton_y/2, ...
        kapton_z_pos - kapton_z/2 : kapton_z_pos + kapton_z/2) = 1;

    % define kapton properties
    KAPTON_DENSITY = 1430;          % [kg/m^3]
    SOUND_SPEED_KAPTON = 2730;      % [m/s]
    
    medium.density(Kapton_foil == 1) =  KAPTON_DENSITY;
    medium.sound_speed(Kapton_foil == 1) = SOUND_SPEED_KAPTON;

    % define air boundary
    air_x_pos = Nx/2;
    air_y_pos = Ny/2;
    air_z_pos = 5;
    %air_z_pos = ((Nz/2) - kapton_z_pos)/2 +2;
    Air = zeros(Nx, Ny, Nz);

    air_x = 150;           % [grid points], 15 mm
    air_y = 150;           % [grid points], 15 mm
    air_z = 9;             % [grid points], 10 mm

    Air(air_x_pos - air_x/2 : air_x_pos + air_x/2, ...
        air_y_pos - air_y/2 : air_y_pos+ air_y/2, ...
        air_z_pos - air_z/2 : air_z_pos + air_z/2) = 1;

    AIR_DENSITY = 1.204;            % [kg/m^3]
    SOUND_SPEED_AIR = 343;          % [m/s]

    medium.density(Air == 1) =  AIR_DENSITY;
    medium.sound_speed(Air == 1) = SOUND_SPEED_AIR;

end


%% =========================================================================
% DEFINE SENSOR ARRAY
% =========================================================================

% create empty kWaveArray object to hold the sensor elements
my_sensor = kWaveArray;

if array_transducer

    %  define some random array element positions 
    N_elements = NUMBER_OF_ELEMENTS;                    
    array_radius = ARRAY_RADIUS;                        % array radius [m]
    array_z_offset = ARRAY_Z_OFFSET;                    % array z-offset [m]
    element_diameter = ARRAY_ELEMENT_DIAMETER;          % element diamater [m]
   
    %%%% golden spiral - evenly spaced elements using k-Wave function %%%%
    % N_elements*2 makes sure all elements are places on hemisphere and not sphere
    sphere = makeCartSphere(array_radius, N_elements*2, [0,0,0], true);
    
    %%% remove negative z-values to create hemisphere %%%
    x = sphere(1,:);
    y = sphere(2,:);
    z = sphere(3,:);
    
    idx = find(z < 0); % finds index of negative z-values
    
    % removes elements with negative z-values
    x( :, idx) = [];
    y( :, idx) = [];
    z( :, idx) = [];
    % merge arrays back to a 3D matrix

    elements = reshape([x y z-array_z_offset].',[],3);
    elements = elements.';
    sensor_positions = elements;
    
    %figure;
    %plot3(elements(1,:),elements(2,:),elements(3,:));
    figure;
    plot3(elements(1,:),elements(2,:),elements(3,:), '.', 'MarkerSize', 50);
    xlim([-0.01 0.01])
    ylim([-0.01 0.01])
    zlim([-0.002 0.008])
    xlabel('[m]');
    ylabel('[m]');
    zlabel('[m]');
    view(90,0)
    f = gcf;
    exportgraphics(f, ['sensor_array_side', '.png'], 'Resolution', 300)
    
    figure;
    plot3(elements(1,:),elements(2,:),elements(3,:), '.', 'MarkerSize', 50);
    xlim([-0.01 0.01])
    ylim([-0.01 0.01])
    zlim([-0.001 0.01])
    xlabel('[m]');
    ylabel('[m]');
    zlabel('[m]');
    view(0,90)
    f = gcf;
    exportgraphics(f, ['sensor_array_top', '.png'], 'Resolution', 300)

    N_elements = length(sensor_positions);
    
    % Remove sensor positions separated by < element_diameter, or by
    % <sqrt(2)*dx so that we have each detector is separated by at least a grid point. 
    closest_distance = 1.01*max([element_diameter, sqrt(2)*kgrid.dx]);
    distances_between_elements = squareform(pdist(sensor_positions'));
    [element_index,dummy] = find((distances_between_elements + eye(N_elements)) < closest_distance);
    sensor_positions(:,element_index')=[];
    
    % the current number of suitably separated elements
    N_elements_separated = size(sensor_positions,2);
    
    % Define the array elements in the k-Wave array using the off-grid-sources
    % approach
    for loop_array = 1:N_elements_separated  
        normal_towards = [0, 0, array_z_offset];
        my_sensor.addDiscElement(sensor_positions(:,loop_array), element_diameter, normal_towards);
    end

end

if single_element_transducer
    
    %element_diameter = SINGLE_ELEMENT_DIAMETER;
    element_diameter = ARRAY_ELEMENT_DIAMETER;
    sensor_positions = [0; 0; 0.006];

    % Define the array elements in the k-Wave array using the off-grid-sources approach
    normal_towards = [0, 0, 0.001];
    my_sensor.addDiscElement(sensor_positions, element_diameter, normal_towards);
end

if matrix_array 
    % Define the pitch
    pitch = 0.3e-3;  % 0.3 mm

    % Number of elements in each bank
    N_banks = 4;
    elements_per_bank_x = 8;
    elements_per_bank_y = 32;
    
    % Focus position 
    focus_pos = [0, 0, 0];
    
    % Initialize an empty matrix to store sensor positions
    sensor_positions = [];

    % Loop through each bank
    for bank_x = 1:elements_per_bank_x
        for bank_y = 1:elements_per_bank_y
            for bank = 1:N_banks
    
                % to center the array in the x-y plane
                x_offset = 0.0051;
                y_offset = 0.0047;

                % Calculate the element center position in the grid
                element_x = (bank_x - 1) * pitch + (bank - 1) * 9 * pitch - x_offset;
                element_y = (bank_y - 1) * pitch - y_offset;

                z_offset = 0.005;
                element_z = z_offset;    % All elements are at the same z-coordinate

                % Append the current sensor position to the matrix
                if on_Axis
                    sensor_positions = [sensor_positions; element_x, element_y, element_z];
                end
                if Side
                    sensor_positions = [sensor_positions; -element_z, element_x, element_y];
                end
              
                % Add square elements to the sensor
                Lx = 275*1e-6;
                Ly = Lx;
                theta = [0,0,0];
                if on_Axis
                    my_sensor.addRectElement([element_x, element_y, element_z], Lx, Ly, theta);
                end
                if Side
                    my_sensor.addRectElement([-element_z, element_x, element_y], Lx, Ly, theta);
                end

            end
        end
    end 
    sensor_positions = sensor_positions.';
end


if column_array
    % Define the pitch
    pitch = 0.2e-3;  % 0.2 mm

    % Number of elements in each bank
    N_banks = 1;
    elements_per_bank_x = 128;
    elements_per_bank_y = 2;
    
    % Focus position 
    focus_pos = [0, 0, 0];
    
    % Initialize an empty matrix to store sensor positions
    sensor_positions = [];

    % Loop through each bank
    for bank_x = 1:elements_per_bank_x
        for bank_y = 1:elements_per_bank_y
            for bank = 1:N_banks
    
                % to center the array in the x-y plane
                x_offset = 0.013;
                y_offset = 0;

                % Calculate the element center position in the grid
                element_x = (bank_x - 1) * pitch + (bank - 1) * 9 * pitch - x_offset;
                element_y = (bank_y - 1) * pitch - y_offset;

                z_offset = 0.005;
                element_z = z_offset;    % All elements are at the same z-coordinate

                % Append the current sensor position to the matrix
                sensor_positions = [sensor_positions; element_x, element_y, element_z];
    
                % Add square elements to the sensor
                Lx = 0.175*1e-3;
                Ly = 25.6*1e-3;
                theta = [0,0,0];
                my_sensor.addRectElement([element_x, element_y, element_z], Lx, Ly, theta);
            end
        end
    end 
    sensor_positions = sensor_positions.';
end


if linear_array
    % Define the pitch
    pitch = 0.23e-3;  % 0.2 mm

    % Number of elements in each bank
    N_banks = 1;
    elements_per_bank_x = 192;
    elements_per_bank_y = 1;
    
    % Focus position 
    focus_pos = [0, 0, 0];
    
    % Initialize an empty matrix to store sensor positions
    sensor_positions = [];

    % Loop through each bank
    for bank_x = 1:elements_per_bank_x
        for bank_y = 1:elements_per_bank_y
            for bank = 1:N_banks
    
                % to center the array in the x-y plane
                x_offset = 0.022;
                y_offset = 0;

                % Calculate the element center position in the grid
                element_x = (bank_x - 1) * pitch + (bank - 1) * 9 * pitch - x_offset;
                element_y = (bank_y - 1) * pitch - y_offset;

                z_offset = 0.005;
                element_z = z_offset;    % All elements are at the same z-coordinate

                % Append the current sensor position to the matrix
                sensor_positions = [sensor_positions; element_x, element_y, element_z];
    
                % Add square elements to the sensor
                Lx = 205*1e-6;
                Ly = Lx;
                theta = [0,0,0];
                my_sensor.addRectElement([element_x, element_y, element_z], Lx, Ly, theta);
            end
        end
    end 
    sensor_positions = sensor_positions.';
end

%figure;
% Plot the array
%my_sensor.plotArray
%grid on
%f = gcf;
%exportgraphics(f, ['sensor_geometry', '.png'], 'Resolution', 300);

figure;
% Plot the front face of the sensor
if matrix_array
    plot(sensor_positions(1, :)*1e3, sensor_positions(2, :)*1e3, 's', ...
        'MarkerSize', 9, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor','blue');
    xlabel('x [mm]');
    ylabel('y [mm]');
    xlim([-6, 6]); 
    ylim([-6, 6]); 
    axis equal;
    f = gcf;
    exportgraphics(f, ['matrix_sensor_geometry', '.png'], 'Resolution', 300);
end

if column_array
    plot(sensor_positions(1, :)*1e3, sensor_positions(2, :)*1e3, 's', ...
        'MarkerSize', 2, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor','blue');
    xlabel('x [mm]');
    ylabel('y [mm]');
    %xlim([-6, 6]); 
    %ylim([-6, 6]);
    axis equal;
    f = gcf;
    exportgraphics(f, ['column_sensor_geometry', '.png'], 'Resolution', 300);
end

if linear_array
    plot(sensor_positions(1, :)*1e3, sensor_positions(2, :)*1e3, 's', ...
        'MarkerSize', 2, 'MarkerFaceColor', 'blue', 'MarkerEdgeColor','blue');
    xlabel('x [mm]');
    ylabel('y [mm]');
    %xlim([-6, 6]); 
    %ylim([-6, 6]);
    axis equal;
    f = gcf;
    exportgraphics(f, ['linear_sensor_geometry', '.png'], 'Resolution', 300);
end


% Simple visual check that the array is aligned around the Bragg peak
% (ignore the odd patterning around the array elements - this is a
% consequence of the off-grid-sources approach)
array_weights = my_sensor.getArrayGridWeights(kgrid);

figure;
%sgtitle('Experiment');

subplot(2,2,1);
%imagesc(squeeze(sum(array_weights,3)./max(max(sum(array_weights,3))) + sum(source.p0,3)./max(max(sum(source.p0,3))) ))
imagesc(y_vec*1e3,x_vec*1e3,squeeze(sum(array_weights,3)./max(max(sum(array_weights,3))) + sum(source.p0,3)./max(max(sum(source.p0,3))) ))
title('x-y plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');

subplot(2,2,2);
%imagesc(squeeze(sum(array_weights,2)./max(max(sum(array_weights,2))) + sum(source.p0,2)./max(max(sum(source.p0,2))) ))
imagesc(y_vec*1e3,x_vec*1e3,squeeze(sum(array_weights,2)./max(max(sum(array_weights,2))) + sum(source.p0,2)./max(max(sum(source.p0,2))) ))
title('x-z plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');

subplot(2,2,3);
%imagesc(squeeze(sum(array_weights)./max(max(sum(array_weights))) + sum(source.p0)./max(max(sum(source.p0))) ))
imagesc(y_vec*1e3,x_vec*1e3,squeeze(sum(array_weights)./max(max(sum(array_weights))) + sum(source.p0)./max(max(sum(source.p0))) ))
title('y-z plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');

f = gcf;
if array_transducer
    exportgraphics(f, ['golden_spiral_array', '.png'], 'Resolution', 300)
end

if single_element_transducer
    exportgraphics(f, ['single_element', '.png'], 'Resolution', 300)
end

if matrix_array
    exportgraphics(f, ['matrix_array', '.png'], 'Resolution', 300)
end
if column_array
    exportgraphics(f, ['column_array', '.png'], 'Resolution', 300)
end
if linear_array
    exportgraphics(f, ['linear_array', '.png'], 'Resolution', 300)
end


%% ========================================================================
% DEFINE A SCINTILLATING FIBRE PLANE STATIONS
% =========================================================================

% POSITIONS NOT FIXED SINCE WE ARE NOT USING THEM
if SciFiStations

    % define the station location
    station_x_pos = Nx/2;
    station_y_pos = Ny/2;
    
    % define station dimensions
    station_thickness = 2.5;        % [grid points], 0.25 mm
    %station_thickness = 5;         % [grid points], 0.5 mm (2 planes)
    station_x = 179;                % [grid points], 20 mm
    station_y = 179;                % [grid points], 20 mm
    station_z = station_thickness;  % [grid points], 0.5 mm (2 planes)
     
    % define material properties
    FIBRE_DENSITY = 1050;       % [kg/m^3]
    BULK_MODULUS = 3250*1000;   % [Pa]
    SOUND_SPEED_FIBRE = sqrt(BULK_MODULUS/FIBRE_DENSITY);  % [m/s]

    % Station 1
    SciFi_Station1 = zeros(Nx, Ny, Nz);
    station1_z_pos = Nz/2 - 30;
    SciFi_Station1(station_x_pos - station_x/2 : station_x_pos + station_x/2, ...
        station_y_pos - station_y/2 : station_y_pos+ station_y/2, ...
        station1_z_pos - station_z/2 : station1_z_pos + station_z/2) = 1;
    medium.density(SciFi_Station1 == 1) = FIBRE_DENSITY;
    medium.sound_speed(SciFi_Station1 == 1) = SOUND_SPEED_FIBRE;

    % Station 2
    SciFi_Station2 = zeros(Nx, Ny, Nz);
    station2_z_pos = Nz/2 - 20;
    SciFi_Station2(station_x_pos - station_x/2 : station_x_pos + station_x/2, ...
        station_y_pos - station_y/2 : station_y_pos+ station_y/2, ...
        station2_z_pos - station_z/2 : station2_z_pos + station_z/2) = 1;
    medium.density(SciFi_Station2 == 1) = FIBRE_DENSITY;
    medium.sound_speed(SciFi_Station2 == 1) = SOUND_SPEED_FIBRE;

    % Station 3
    SciFi_Station3 = zeros(Nx, Ny, Nz);
    station3_z_pos = Nz/2 - 10;
    SciFi_Station3(station_x_pos - station_x/2 : station_x_pos + station_x/2, ...
        station_y_pos - station_y/2 : station_y_pos+ station_y/2, ...
        station3_z_pos - station_z/2 : station3_z_pos + station_z/2) = 1;
    medium.density(SciFi_Station3 == 1) = FIBRE_DENSITY;
    medium.sound_speed(SciFi_Station3 == 1) = SOUND_SPEED_FIBRE;

    % Station 4 
    SciFi_Station4 = zeros(Nx, Ny, Nz);
    station4_z_pos = Nz/2;
    SciFi_Station4(station_x_pos - station_x/2 : station_x_pos + station_x/2, ...
        station_y_pos - station_y/2 : station_y_pos+ station_y/2, ...
        station4_z_pos - station_z/2 : station4_z_pos + station_z/2) = 1;
    medium.density(SciFi_Station4 == 1) = FIBRE_DENSITY;
    medium.sound_speed(SciFi_Station4 == 1) = SOUND_SPEED_FIBRE;

end


%% ========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the time limit to belong enough to capture the wave
t_end = sqrt(2)*Nx*dx / SOUND_SPEED_WATER;   % [s]
%t_end = t_end*2;   % to capture the reflected wave too (if there is)

% set the timestep (smaller CFL -> smaller timestep)
kgrid.makeTime(medium.sound_speed, CFL, t_end);
filename = append('kgrid_t_array_',add,'.csv');
writematrix(kgrid.t_array, filename); 


% assign binary mask in order to be able to run the simulation
% (in future this step will be integrated into k-Wave, and we will be able
% to simply pass my_sensor)
sensor.mask = my_sensor.getArrayBinaryMask(kgrid);

if Kapton_window
    if SciFiStations
        display_mask = Kapton_foil + SciFi_Station1 + SciFi_Station2 + ...
            SciFi_Station3 + SciFi_Station4;
    else
        display_mask = Kapton_foil;
        %display_mask = Kapton_foil + Air;
        %display_mask = Kapton_foil + Air + sensor.mask;
    end
else
    if SciFiStations
        display_mask = SciFi_Station1 + SciFi_Station2 + ...
            SciFi_Station3 + SciFi_Station4;
    else
        display_mask = 'off';
    end
end

movie_name = append('wave_movie_',add);

% absorption coefficient for the PML, 0 = full reflective, 1 = full absorbtive
PML_alpha = 1-0.84;  % i.e., 84% reflective

input_args = {'PMLInside', false, 'PMLSize',PML_SIZE, 'PMLAlpha', PML_alpha, 'PlotPML', ...
        false, 'PlotSim', false, 'DataCast', 'single', 'Smooth', false, 'PlotScale', ...
        [-max(source.p0(:)) max(source.p0(:))]/2, 'RecordMovie', false, ...
        'MovieName', movie_name, 'DisplayMask', display_mask};
    
% run the simulation
%sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});    % runs in native matlab
sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});  % uses the compiled C++ code
%sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});  % used the compiled CUDA code

% apvvplying sensor response
Fs = 1/kgrid.dt;  % [Hz]
if SensorResponse
    sensor_data = gaussianFilter(sensor_data, Fs, center_freq, bandwidth);
else 
    sensor_data = sensor_data;
end

% Add electronic noise to each sensor element with a different seed
if Electronic_noise
    % Use 'shuffle' to set the seed based on the current time
    rng('shuffle');
    st_dev = 0.1;

    % Generate a shuffled array of seeds for each sensor element
    num_elements = numel(sensor_data);
    shuffled_indices = randperm(num_elements);
    seeds = shuffled_indices;  % Using shuffled indices as seeds

    % Initialize an array to store noisy sensor data
    noisy_sensor_data = zeros(size(sensor_data));

    % Introduce random noise to each sensor element using a different seed
    for i = 1:size(sensor_data, 1)
        for j = 1:size(sensor_data, 2)
            rng(seeds((i-1)*size(sensor_data, 2) + j));  % Set the seed for each element
            noisy_sensor_data(i, j) = sensor_data(i, j) + st_dev * randn();  
        end
    end
    sensor_data = noisy_sensor_data;
end


if array_transducer
    % combine data to give one time series per array element
    % (in future this step will be integrated into k-Wave)
    combined_sensor_data = my_sensor.combineSensorData(kgrid, sensor_data);

    % sums pressures at every detector voxel for each time step 
    combined = sum(combined_sensor_data, 1); 
    
    [fr, sensor_data_as] = spect(combined, 1/kgrid.dt, 'Dim', 2, 'Plot', false);
    
    filename = append('sensor_data_sum_array_',add,'.csv');
    writematrix(combined, filename);
    filename = append('amplitude_data_array_',add,'.csv');
    writematrix(sensor_data_as, filename);
    filename = append('freq_data_array_',add,'.csv');          % [Hz]
    writematrix(fr, filename); 
end


if single_element_transducer   

    % sums pressures at every detector voxel for each time step 
    combined = sum(sensor_data, 1); 
    [fr, sensor_data_as] = spect(combined, 1/kgrid.dt, 'Dim', 2, 'Plot', false);
    %[fr, sensor_data_as] = spect(sensor_data, 1/kgrid.dt, 'Dim', 2, 'Plot', false);
    
    filename = append('sensor_data_sum_single_element_',add,'.csv');
    writematrix(combined, filename);
    filename = append('amplitude_data_single_element_',add,'.csv');
    writematrix(sensor_data_as, filename);
    filename = append('freq_data_single_element_',add,'.csv');  % [Hz]
    writematrix(fr, filename); 
end

if matrix_array
    % combine data to give one time series per array element
    % (in future this step will be integrated into k-Wave)
    combined_sensor_data = my_sensor.combineSensorData(kgrid, sensor_data);

    % sums pressures at every detector voxel for each time step 
    combined = sum(combined_sensor_data, 1); 
    
    [fr, sensor_data_as] = spect(combined, 1/kgrid.dt, 'Dim', 2, 'Plot', false);
    
    filename = append('sensor_data_sum_matrix_array_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(combined, filename);
    filename = append('amplitude_data_matrix_array_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(sensor_data_as, filename);
    filename = append('freq_data_matrix_array_',add,add_Response,add_ElectronicNoise,'.csv');      % [Hz]
    writematrix(fr, filename); 
end

if column_array
    % combine data to give one time series per array element
    % (in future this step will be integrated into k-Wave)
    combined_sensor_data = my_sensor.combineSensorData(kgrid, sensor_data);

    % sums pressures at every detector voxel for each time step 
    combined = sum(combined_sensor_data, 1); 
    
    [fr, sensor_data_as] = spect(combined, 1/kgrid.dt, 'Dim', 2, 'Plot', false);
    
    filename = append('sensor_data_sum_column_array_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(combined, filename);
    filename = append('amplitude_data_column_array_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(sensor_data_as, filename);
    filename = append('freq_data_column_array_',add,add_Response,add_ElectronicNoise,'.csv');      % [Hz]
    writematrix(fr, filename); 
end

if linear_array
    % combine data to give one time series per array element
    % (in future this step will be integrated into k-Wave)
    combined_sensor_data = my_sensor.combineSensorData(kgrid, sensor_data);

    % sums pressures at every detector voxel for each time step 
    combined = sum(combined_sensor_data, 1); 
    
    [fr, sensor_data_as] = spect(combined, 1/kgrid.dt, 'Dim', 2, 'Plot', false);
    
    filename = append('sensor_data_sum_linear_array_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(combined, filename);
    filename = append('amplitude_data_linear_array_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(sensor_data_as, filename);
    filename = append('freq_data_linear_array_',add,add_Response,add_ElectronicNoise,'.csv');      % [Hz]
    writematrix(fr, filename); 
end


%% ========================================================================
% IMAGE RECONSTRUCTION (TIME REVERSAL)
% =========================================================================

% removes the initial pressure field used in the simulation
source = rmfield(source, 'p0');

% uses the sensor points as sources in time reversal
% (associate each time series to the detector's centre point as a first go)
source.p_mask = cart2grid(kgrid, sensor_positions);

if array_transducer
    % times reverse and assign the data
    source.p = fliplr(combined_sensor_data);
end
if single_element_transducer
    % times reverse and assign the data
    source.p = fliplr(sensor_data);
end
if matrix_array
    % times reverse and assign the data
    source.p = fliplr(combined_sensor_data);
end
if column_array
    % times reverse and assign the data
    source.p = fliplr(combined_sensor_data);
end
if linear_array
    % times reverse and assign the data
    source.p = fliplr(combined_sensor_data);
end

% enforce, rather than add, the time-reversed pressure values
source.p_mode = 'dirichlet';  

% sets the simulation to record the final image (at t = 0)
sensor.record = {'p_final'};

% runs the time reversal reconstruction
%p0_estimate = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
p0_estimate = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
%p0_estimate = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});

% applies a positivity condition
% factor of 2 applied for first order correction due to use of hemisphere
p0_estimate.p_final = p0_estimate.p_final .* (p0_estimate.p_final > 0);

% stores the latest image estimate
p0_1 = p0_estimate.p_final;


%% =========================================================================
% ITERATE THE TR RECONSTRUCTION TO IMPROVE THE IMAGE
% =========================================================================

for loop = 1:NUMBER_OF_TR_ITERATIONS
    
    % removes the source used in the previous time reversal
    source = rmfield(source, 'p');

    % sets the initial pressure to be the latest estimate of p0
    source.p0 = p0_estimate.p_final;
    
    % set sthe simulation to record the time series
    sensor = rmfield(sensor, 'record');
    
    % calculates the time series using the latest estimate of p0
    %sensor_data2 = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
    sensor_data2 = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:}); 
    %sensor_data2 = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
    
    %combined_sensor_data = my_sensor.combineSensorData(kgrid, sensor_data2);
    % sums pressures at every detector voxel for each time step 
    %combined = sum(combined_sensor_data, 1); 

    % calculates the error in the estimated time series
    data_difference = sensor_data - sensor_data2;
    
    % assigns the data_difference as a time-reversal source
    source.p_mask = sensor.mask;
    source.p = fliplr(data_difference);
    source = rmfield(source, 'p0');
    source.p_mode = 'dirichlet';   

    % sets the simulation to record the final image (at t = 0)
    sensor.record = {'p_final'};
    
    % runs the time reversal reconstruction
    %p0_update = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
    p0_update = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
    %p0_update = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
    
    % adds the update to the latest image    
    p0_estimate.p_final = p0_estimate.p_final + p0_update.p_final;

    % applies a positivity condition
    % factor of 2 applied for first order correction due to use of hemisphere
    p0_estimate.p_final = p0_estimate.p_final .* (p0_estimate.p_final > 0); 
    
    % stores the latest image estimate
    eval(['p0_' num2str(loop) ' = p0_estimate.p_final;']);

end

%% =========================================================================
% VISUALISATION
% =========================================================================

%padding_added = [20,20,20];
%p0_1 = padarray(p0_1, padding_added);
%p0_2 = padarray(p0_2, padding_added);
%p0_3 = padarray(p0_3, padding_added);

%Nx = 200;
%Ny = 200;
%Nz = 200;

%%%% integrating radially the reconstruction %%%%

%%%% R1 %%%%
% empty array to store values
int_x_p = zeros(1,Nx);
int_y_p = zeros(1,Ny);
int_z_p = zeros(1,Nz);

plane_xy_int = zeros(Nz);
plane_xz_int = zeros(Ny);
plane_yz_int = zeros(Nx);



% integrate along each plane
% and store total pressure in each plane - along z-axi
for p = 1:Nx
    
    % pressure in each plane = 1 voxel
    plane_xy = squeeze(p0_1(:,:,p));
    plane_xz = squeeze(p0_1(:,p,:));
    plane_yz = squeeze(p0_1(p,:,:));
    
    % adding the pressure in all voxels
    combined_xy = sum(plane_xy, 'all'); 
    combined_xz = sum(plane_xz, 'all'); 
    combined_yz = sum(plane_yz, 'all'); 
    % store into arrays
    int_x_p(p) = combined_yz;
    int_y_p(p) = combined_xz;
    int_z_p(p) = combined_xy;

    % adds the images
    plane_xy_int = plane_xy_int + plane_xy;
    plane_xz_int = plane_xz_int + plane_xz;
    plane_yz_int = plane_yz_int + plane_yz;
end

if array_transducer
    filename = append('array_int_p0_z_rec_1_',add,add_Response,'.csv');
    writematrix(int_z_p, filename); 
end
if single_element_transducer
    filename = append('single_element_int_p0_z_rec_1_',add,add_Response,'.csv');
    writematrix(int_z_p, filename); 
end
if matrix_array
    filename = append('matrix_array_int_p0_x_rec_1_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_x_p, filename);
    filename = append('matrix_array_int_p0_y_rec_1_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_y_p, filename);
    filename = append('matrix_array_int_p0_z_rec_1_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_z_p, filename); 
end
if column_array
    filename = append('column_array_int_p0_z_rec_1_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_z_p, filename); 
end
if linear_array
    filename = append('linear_array_int_p0_z_rec_1_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_z_p, filename); 
end

%%%% R2 %%%%
% empty array to store values
int_x_p = zeros(1,Nx);
int_y_p = zeros(1,Ny);
int_z_p = zeros(1,Nz);

plane_xy_int = zeros(Nz);
plane_xz_int = zeros(Ny);
plane_yz_int = zeros(Nx);

% integrate along each plane
% and store total pressure in each plane - along z-axi
for p = 1:Nx
    
    % pressure in each plane = 1 voxel
    plane_xy = squeeze(p0_2(:,:,p));
    plane_xz = squeeze(p0_2(:,p,:));
    plane_yz = squeeze(p0_2(p,:,:));
    
    % adding the pressure in all voxels
    combined_xy = sum(plane_xy, 'all'); 
    combined_xz = sum(plane_xz, 'all'); 
    combined_yz = sum(plane_yz, 'all'); 
    % store into arrays
    int_x_p(p) = combined_yz;
    int_y_p(p) = combined_xz;
    int_z_p(p) = combined_xy;

    % adds the images
    plane_xy_int = plane_xy_int + plane_xy;
    plane_xz_int = plane_xz_int + plane_xz;
    plane_yz_int = plane_yz_int + plane_yz;
end

if array_transducer
    filename = append('array_int_p0_z_rec_2_',add,add_Response,'.csv');
    writematrix(int_z_p, filename); 
end
if single_element_transducer
    filename = append('single_element_int_p0_z_rec_2_',add,add_Response,'.csv');
    writematrix(int_z_p, filename); 
end
if matrix_array
    filename = append('matrix_array_int_p0_x_rec_2_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_x_p, filename); 
    filename = append('matrix_array_int_p0_y_rec_2_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_y_p, filename); 
    filename = append('matrix_array_int_p0_z_rec_2_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_z_p, filename); 
end
if column_array
    filename = append('column_array_int_p0_z_rec_2_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_z_p, filename); 
end
if linear_array
    filename = append('linear_array_int_p0_z_rec_2_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_z_p, filename); 
end
 

%%%% R3 %%%%
% empty array to store values
int_x_p = zeros(1,Nx);
int_y_p = zeros(1,Ny);
int_z_p = zeros(1,Nz);

plane_xy_int = zeros(Nz);
plane_xz_int = zeros(Ny);
plane_yz_int = zeros(Nx);

% integrate along each plane
% and store total pressure in each plane - along z-axi
for p = 1:Nx
    
    % pressure in each plane = 1 voxel
    plane_xy = squeeze(p0_3(:,:,p));
    plane_xz = squeeze(p0_3(:,p,:));
    plane_yz = squeeze(p0_3(p,:,:));
    
    % adding the pressure in all voxels
    combined_xy = sum(plane_xy, 'all'); 
    combined_xz = sum(plane_xz, 'all'); 
    combined_yz = sum(plane_yz, 'all'); 
    % store into arrays
    int_x_p(p) = combined_yz;
    int_y_p(p) = combined_xz;
    int_z_p(p) = combined_xy;

    % adds the images
    plane_xy_int = plane_xy_int + plane_xy;
    plane_xz_int = plane_xz_int + plane_xz;
    plane_yz_int = plane_yz_int + plane_yz;
end

if array_transducer
    filename = append('array_int_p0_z_rec_3_',add,add_Response,'.csv');
    writematrix(int_z_p, filename); 
end
if single_element_transducer
    filename = append('single_element_int_p0_z_rec_3_',add,add_Response,'.csv');
    writematrix(int_z_p, filename); 
end
if matrix_array
    filename = append('matrix_array_int_p0_x_rec_3_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_x_p, filename); 
    filename = append('matrix_array_int_p0_y_rec_3_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_y_p, filename); 
    filename = append('matrix_array_int_p0_z_rec_3_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_z_p, filename); 
end
if column_array
    filename = append('column_array_int_p0_z_rec_3_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_z_p, filename); 
end
if linear_array
    filename = append('linear_array_int_p0_z_rec_3_',add,add_Response,add_ElectronicNoise,'.csv');
    writematrix(int_z_p, filename); 
end


% find the maximum pressure - comment out to use maximum at source
%maxValue = max([max(plane_xy_int), max(plane_xz_int), max(plane_yz_int)]);

figure;
%sgtitle('Integrated Reconstructed Pressure Distribution')
subplot(2,2,1);
%imagesc(plane_xy_int);
imagesc(y_vec*1e3,x_vec*1e3,plane_xy_int); % plot in [mm]
title('x-y plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
caxis([0, maxValue]);
c.Label.String = 'Pressure [Pa]';

subplot(2,2,2);
%imagesc(plane_xz_int);
imagesc(y_vec*1e3,x_vec*1e3,plane_xz_int); % plot in [mm]
title('x-z plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
caxis([0, maxValue]);
c.Label.String = 'Pressure [Pa]';

subplot(2,2,3);
%imagesc(plane_yz_int);
imagesc(y_vec*1e3,x_vec*1e3,plane_yz_int); % plot in [mm]
title('y-z plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
caxis([0, maxValue]);
c.Label.String = 'Pressure [Pa]';
f = gcf;

if array_transducer
    filename = append('array_integrated_pressure_distribution_reconstructed_',add,add_Response,add_ElectronicNoise);
    exportgraphics(f, [filename, '.png'], 'Resolution', 300);
end
if single_element_transducer
    filename = append('single_element_integrated_pressure_distribution_reconstructed_',add,add_Response,add_ElectronicNoise);
    exportgraphics(f, [filename, '.png'], 'Resolution', 300);
end
if matrix_array
    filename = append('matrix_array_integrated_pressure_distribution_reconstructed_',add,add_Response,add_ElectronicNoise);
    exportgraphics(f, [filename, '.png'], 'Resolution', 300);
end
if column_array
    filename = append('column_array_integrated_pressure_distribution_reconstructed_',add,add_Response,add_ElectronicNoise);
    exportgraphics(f, [filename, '.png'], 'Resolution', 300);
end
if linear_array
    filename = append('linear_array_integrated_pressure_distribution_reconstructed_',add,add_Response,add_ElectronicNoise);
    exportgraphics(f, [filename, '.png'], 'Resolution', 300);
end


% Plot difference image
figure;
%sgtitle('Integrated Reconstructed Pressure Distribution')
subplot(2,2,1);
%imagesc(plane_xy_int-plane_xy_int_source);
imagesc(y_vec*1e3,x_vec*1e3,plane_xy_int-plane_xy_int_source); % plot in [mm]
title('x-y plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
%caxis([0, maxValue]);
c.Label.String = 'Pressure [Pa]';

subplot(2,2,2);
%imagesc(plane_xz_int-plane_xz_int_source);
imagesc(y_vec*1e3,x_vec*1e3,plane_xz_int-plane_xz_int_source); % plot in [mm]
title('x-z plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
%caxis([0, maxValue]);
c.Label.String = 'Pressure [Pa]';

subplot(2,2,3);
%imagesc(plane_yz_int-plane_yz_int_source);
imagesc(y_vec*1e3,x_vec*1e3,plane_yz_int-plane_yz_int_source); % plot in [mm]
title('y-z plane');
%xlabel('[voxels]');
%ylabel('[voxels]');
xlabel('[mm]');
ylabel('[mm]');
colormap(getColorMap);
c = colorbar;
%caxis([0, maxValue]);
c.Label.String = 'Pressure [Pa]';
f = gcf;

filename = append('integrated_pressure_distribution_difference_',add,add_Response,add_ElectronicNoise);
exportgraphics(f, [filename, '.png'], 'Resolution', 300);

toc
    
