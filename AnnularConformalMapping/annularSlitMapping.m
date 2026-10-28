function[tau] = annularSlitMapping(v, f, inner_bd_id)
% Author: Shen Changqing 2024.08.15
% Email: 625372005@qq.com
% Feel free to contact me if you have any questions.
% If you use this code in your own work, please cite the following paper:
% [1] X. Ban, M. Goswami, W. Zeng, X. Gu and J. Gao, "Topology dependent space filling curves for sensor networks and applications," 2013 Proceedings IEEE INFOCOM, Turin, Italy, 2013, pp. 2166-2174, doi: 10.1109/INFCOM.2013.6567019.
% [2] G. P. T. Choi, Efficient Conformal Parameterization of Multiply-Connected Surfaces Using Quasi-Conformal Theory. Journal of Scientific Computing, 87(3), 70, 2021.
% [3] Wang, Y., Gu, X., Chan, T.F., Thompson, P.M., Yau, ST. (2008). Conformal Slit Mapping and Its Applications to Brain Surface Parameterization. In: Metaxas, D., Axel, L., Fichtinger, G., Székely, G. (eds) Medical Image Computing and Computer-Assisted Intervention – MICCAI 2008. MICCAI 2008. Lecture Notes in Computer Science, vol 5241. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-85988-8_70
% [4] Yin, X., Dai, J., Yau, ST., Gu, X. (2008). Slit Map: Conformal Parameterization for Multiply Connected Surfaces. In: Chen, F., Jüttler, B. (eds) Advances in Geometric Modeling and Processing. GMP 2008. Lecture Notes in Computer Science, vol 4975. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-79246-8_31

% rely on:
% insertrows()

[all_v_sliced, all_f_sliced, all_slice_path, outer_bd_id] = sliceMesh(v, f);%有bug

all_f_0 = {};
all_eta = {};
for idx = 1:length(all_slice_path)
    if idx == outer_bd_id
        all_f_0{idx} = 0;
        all_eta{idx} = 0;
        continue;
    end
    all_f_0{idx} = exactHarmonicOneFormForSlitMapping(v, f, outer_bd_id, idx);
    all_eta{idx} = closedHarmonicOneFormForSlitMapping(v, f, all_v_sliced{idx}, all_f_sliced{idx}, all_slice_path{idx});
end

all_omega = {};
for idx = 1:length(all_slice_path)
    if idx == outer_bd_id
        all_omega{idx} = f * 0;
        all_omega{idx + length(all_slice_path)} = f * 0;
        continue;
    end
    all_omega{idx} = all_f_0{idx}(f);
    all_omega{idx + length(all_slice_path)} = all_eta{idx}(all_f_sliced{idx});
end

[tau] = calTau(all_omega, outer_bd_id, inner_bd_id, v, f);
