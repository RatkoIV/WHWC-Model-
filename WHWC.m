% Čišćenje radnog prostora
clear;
clc;

% Učitavanje originalne slike u full spectrum
full_spectrum_image = imread('C:\Users\Bato\Desktop\RADOVI 2024\Digital image correlation\Images\full_spectrum_image.jpg');

% Definisanje talasnih dužina i putanje do slika
wavelengths = [446, 550, 649, 657, 710, 721, 750, 765]; % Talasne dužine u nm
num_wavelengths = length(wavelengths);
wavelength_images = cell(num_wavelengths, 1);

% Učitavanje slika za svaku talasnu dužinu
base_path = 'C:\Users\Bato\Desktop\RADOVI 2024\Digital image correlation\Images\nm\';
for k = 1:num_wavelengths
    filename = [base_path, 'wavelength_', num2str(wavelengths(k)), '.jpg'];
    wavelength_images{k} = imread(filename);
end

% Konverzija slika u grayscale ako su u boji
if size(full_spectrum_image, 3) == 3
    full_spectrum_image = rgb2gray(full_spectrum_image);
end

for k = 1:num_wavelengths
    if size(wavelength_images{k}, 3) == 3
        wavelength_images{k} = rgb2gray(wavelength_images{k});
    end
end

% Konverzija slika u tip double za računanje
full_spectrum_image = double(full_spectrum_image);

for k = 1:num_wavelengths
    wavelength_images{k} = double(wavelength_images{k});
end

% Parametri za DIC (veličina prozora i korak)
window_size = 300;  % Smanjen prozor za detaljniju analizu
step_size = 18;     % Smanjen korak za više podataka

% Parametri za WHWC funkciju
sigma = 1;    % Parametar za harmonične težinske faktore

% Veličina slika
[rows, cols] = size(full_spectrum_image);

% Inicijalizacija matrica za rezultate
correlation_matrix = zeros(num_wavelengths, 1);
std_dev_matrix = zeros(num_wavelengths, 1);
min_corr_matrix = zeros(num_wavelengths, 1);
max_corr_matrix = zeros(num_wavelengths, 1);

% Niz za skladištenje svih korelacionih koeficijenata
all_corr_data = cell(num_wavelengths, 1);

% Nizovi za pamćenje broja prozora po slikama
num_windows = zeros(num_wavelengths, 1);

% DIC algoritam za svaku talasnu dužinu
parfor k = 1:num_wavelengths % Koristimo parfor za paralelizaciju
    current_image = wavelength_images{k};
    num_windows(k) = floor((rows-window_size)/step_size + 1) * floor((cols-window_size)/step_size + 1);
    corr_data = zeros(1, num_windows(k)); % Alokacija memorije za sve korelacione koeficijente
    index = 1;
    
    for i = 1:step_size:(rows-window_size)
        for j = 1:step_size:(cols-window_size)
            % Izvlačenje prozora iz full spectrum slike
            window_ref = full_spectrum_image(i:i+window_size-1, j:j+window_size-1);
            % Izvlačenje prozora iz slike na trenutnoj talasnoj dužini
            window_def = current_image(i:i+window_size-1, j:j+window_size-1);
            % Izračunavanje WHWC korelacionog koeficijenta
            corr_coef = whwc(window_ref, window_def, sigma);
            corr_data(index) = corr_coef;
            index = index + 1;
        end
    end
    
    % Statistika za trenutnu talasnu dužinu
    correlation_matrix(k) = mean(corr_data);
    std_dev_matrix(k) = std(corr_data);
    min_corr_matrix(k) = min(corr_data);
    max_corr_matrix(k) = max(corr_data);
    
    all_corr_data{k} = corr_data; % Čuvanje svih podataka za kasniju analizu
end

% Prikaz rezultata - prosečni korelacioni koeficijenti (linijski grafikoni) sa polinomskom interpolacijom
figure;
p_corr = polyfit(wavelengths, correlation_matrix, 3); % Polinom 3. stepena
y_fit_corr = polyval(p_corr, wavelengths);
plot(wavelengths, correlation_matrix, 'o', wavelengths, y_fit_corr, '-r', 'LineWidth', 2);
title('Avg. WHWC Correlation with 3rd Degree Polynomial');
xlabel('Wavelength (nm)');
ylabel('Avg. WHWC Correlation Coefficient');
grid on;

% Prikaz rezultata - standardna devijacija (linijski grafikoni) sa polinomskom interpolacijom
figure;
p_std = polyfit(wavelengths, std_dev_matrix, 3); % Polinom 3. stepena
y_fit_std = polyval(p_std, wavelengths);
plot(wavelengths, std_dev_matrix, 'o', wavelengths, y_fit_std, '-r', 'LineWidth', 2);
title('Std Dev of WHWC Correlation with 3rd Degree Polynomial');
xlabel('Wavelength (nm)');
ylabel('Std Dev of WHWC Correlation Coefficient');
grid on;

% Prikaz rezultata - minimalni korelacioni koeficijenti (linijski grafikoni) sa polinomskom interpolacijom
figure;
p_min = polyfit(wavelengths, min_corr_matrix, 3); % Polinom 3. stepena
y_fit_min = polyval(p_min, wavelengths);
plot(wavelengths, min_corr_matrix, 'o', wavelengths, y_fit_min, '-r', 'LineWidth', 2);
title('Min WHWC Correlation with 3rd Degree Polynomial');
xlabel('Wavelength (nm)');
ylabel('Min WHWC Correlation Coefficient');
grid on;

% Prikaz rezultata - maksimalni korelacioni koeficijenti (linijski grafikoni) sa polinomskom interpolacijom
figure;
p_max = polyfit(wavelengths, max_corr_matrix, 3); % Polinom 3. stepena
y_fit_max = polyval(p_max, wavelengths);
plot(wavelengths, max_corr_matrix, 'o', wavelengths, y_fit_max, '-r', 'LineWidth', 2);
title('Max WHWC Correlation with 3rd Degree Polynomial');
xlabel('Wavelength (nm)');
ylabel('Max WHWC Correlation Coefficient');
grid on;

% Prikaz korelacionih koeficijenata za svaku talasnu dužinu - zasebni linijski grafikoni sa polinomskom interpolacijom
p_ind_polynomials = cell(num_wavelengths, 1); % Čuvanje polinoma za svaku talasnu dužinu
y_fit_ind_polynomials = cell(num_wavelengths, 1); % Čuvanje interpoliranih vrednosti

for k = 1:num_wavelengths
    figure; % Otvara novi prozor za svaku talasnu dužinu
    x = 1:num_windows(k);
    y = all_corr_data{k};
    p_ind = polyfit(x, y, 3); % Polinom 3. stepena
    y_fit_ind = polyval(p_ind, x); % Interpolacija sa polinomom 3. stepena
    p_ind_polynomials{k} = p_ind; % Čuvanje koeficijenata polinoma
    y_fit_ind_polynomials{k} = y_fit_ind; % Čuvanje interpoliranih vrednosti
    plot(x, y, 'o', x, y_fit_ind, '-r', 'LineWidth', 1.5);
    title(['WHWC Correlation for ', num2str(wavelengths(k)), ' nm with 3rd Degree Polynomial']);
    xlabel('Window Index');
    ylabel('WHWC Correlation Coefficient');
    grid on;
end

% Generisanje generalizovane polinomske funkcije za sve talasne dužine
general_x = []; % Svi indeksi prozora
general_y = []; % Svi korelacioni koeficijenti

% Prikupljanje podataka za sve talasne dužine
for k = 1:num_wavelengths
    x = (1:num_windows(k))'; % Indeksi prozora
    y = all_corr_data{k}; % Korelacioni koeficijenti
    
    % Kombinovanje podataka
    general_x = [general_x; x];
    general_y = [general_y; y(:)]; % Uzimanje svih vrednosti y
end

% Provera da li su vektori iste dužine
if length(general_x) == length(general_y)
    % Generisanje polinoma 3. stepena za sve talasne dužine
    p_general = polyfit(general_x, general_y, 3); % Polinom 3. stepena
    y_fit_general = polyval(p_general, general_x);
    
    % Prikaz generalizovanog polinoma
    figure;
    plot(general_x, general_y, 'o', general_x, y_fit_general, '-r', 'LineWidth', 2);
    title('General WHWC Correlation with 3rd Degree Polynomial');
    xlabel('Window Index');
    ylabel('WHWC Correlation Coefficient');
    grid on;

    % Ispis koeficijenata polinoma
    fprintf('Coefficients of the General 3rd-Degree Polynomial:\n');
    disp(p_general);
else
    error('General vectors x and y are not the same length. Please check the data.');
end

% Definicija WHWC funkcije (Windowed Harmonic Weighted Correlation)
function rho_WHWC = whwc(I1, I2, sigma)
    [m, n] = size(I1); % Pretpostavljamo da su I1 i I2 iste veličine

    % Definisanje harmoničnog težinskog faktora W_h(i,j)
    [X, Y] = meshgrid(1:n, 1:m);
    X = double(X); % Konverzija X u tip double
    Y = double(Y); % Konverzija Y u tip double
    W_h = 1 ./ (1 + ((X-n/2).^2 + (Y-m/2).^2) / sigma^2);

    % Sinusni faktor koji uvodi nelinearnost u zavisnosti od indeksa prozora
    Sin_factor = sin(pi/2 * (X .* Y) / (m * n));

    % Izračunavanje WHWC (Windowed Harmonic Weighted Correlation)
    numerator = sum(sum(W_h .* (I1 .* I2) .* Sin_factor));
    denominator = sqrt(sum(sum(W_h .* (I1.^2)))) * sqrt(sum(sum(W_h .* (I2.^2))));

    rho_WHWC = numerator / denominator;
end
