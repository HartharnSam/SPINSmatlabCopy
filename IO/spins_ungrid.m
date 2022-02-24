% unmap
mkdir('remapped')
[x z] = spinsgrid2d;
Nz2 = 256*2;
Nx2 = 2048;
zs = NaN(Nx2, Nz2);
u2 = zs*NaN;
w2 = u2; rho2 = u2;

for i = 0:130
    rho = spins_reader_new('rho', i);
    
    u = spins_reader_new('u', i);
    w = spins_reader_new('w', i);
    for ii = 1:length(x)
        zs(ii, :) = linspace(min(z(ii, :)), max(z(ii, :)), Nz2);
        rho2(ii, :) = interp1(z(ii, :), rho(ii, :), zs(ii, :));
        u2(ii,:) = interp1(z(ii, :), u(ii, :), zs(ii, :));
        w2(ii, :) = interp1(z(ii, :), w(ii, :), zs(ii, :));
    end
    spins_writer(['remapped/rho.', num2str(i)], rho2);
    
    spins_writer(['remapped/u.', num2str(i)], u2);
    spins_writer(['remapped/w.', num2str(i)], w2);
    disp(num2str(i));
end
spins_writer(['remapped/zgrid'], zs);
for jj = 1:Nz2
    xs(:, jj) = x(:, 1);    
end
spins_writer(['remapped/xgrid'], xs);

%copyfile('xgrid', './remapped/xgrid');
copyfile('xgrid_reader.m', './remapped/');
copyfile('zgrid_reader.m', './remapped/');
copyfile('spins.conf', './remapped/');
