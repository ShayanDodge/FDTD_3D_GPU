clc
clear all
close all
uiopen('C:\Users\Mahshad\Desktop\Multilayer\plan_three_layer_non_dispersive_skin_fdtd_cpml_3D_gaussian\three_layer_non_dispersive_skin_fdtd_cpml_gaussian_3D_3D.fig',1)
for i=0:51
    saveas(gcf,['three_layer_non_dispersive_skin_fdtd_cpml_gaussian_3D_3D_' num2str(52-i) '.tiff'])
    close(gcf)
end