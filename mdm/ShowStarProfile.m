%%
i = 20;
mode = 'training';

img = ImagePair.Load(mode, i); 
figure(1);
img.Show();

model = ModelF();
model.ParameterVector = estimate.(mode).X(i, :);

%
nu = img.S.ValidNu;
r = model.PSF.R(nu);
p = model.PSF.ComputeAt(nu);

figure(2);
hold off;
y = img.S.ValidF / model.ScaleS ./ exp(-2i * pi * nu * model.OffsetS.');
plot(r, real(y), '.');
grid on;
hold on;
plot(r, p, 'r.')

%
nu_g = img.G.ValidNu;
r_g = model.Galaxy.R(nu_g);
p_g = model.Galaxy.ComputeAt(nu_g);

figure(3);
hold off;
y_g = img.G.ValidF / model.ScaleG ./ exp(-2i * pi * nu_g * model.OffsetG.');
plot(r_g, real(y_g), '.');
grid on;
hold on;
plot(r_g, p_g .* model.PSF.ComputeAt(nu_g), 'r.');

%

figure(4);
model.GenerateImage(false).Show();