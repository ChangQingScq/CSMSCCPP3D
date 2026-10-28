function[v_sliced, f_sliced, f_0, eta, lambda] = annularConformalMapping(v, f)
% Author: ChangqingScq 2024.08.05
% Email: 625372005@qq.com
%% Modified Based On:
% [1] G. P. T. Choi, Efficient Conformal Parameterization of Multiply-Connected Surfaces Using Quasi-Conformal Theory. Journal of Scientific Computing, 87(3), 70, 2021.
% Feel free to contact me if you have any questions.
% If you use this code in your own work, please cite the following paper:
% [2] X. Ban, M. Goswami, W. Zeng, X. Gu and J. Gao, "Topology dependent space filling curves for sensor networks and applications," 2013 Proceedings IEEE INFOCOM, Turin, Italy, 2013, pp. 2166-2174, doi: 10.1109/INFCOM.2013.6567019.
% [3] Wang, Y., Gu, X., Chan, T.F., Thompson, P.M., Yau, ST. (2008). Conformal Slit Mapping and Its Applications to Brain Surface Parameterization. In: Metaxas, D., Axel, L., Fichtinger, G., Székely, G. (eds) Medical Image Computing and Computer-Assisted Intervention – MICCAI 2008. MICCAI 2008. Lecture Notes in Computer Science, vol 5241. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-85988-8_70
% [4] Yin, X., Dai, J., Yau, ST., Gu, X. (2008). Slit Map: Conformal Parameterization for Multiply Connected Surfaces. In: Chen, F., Jüttler, B. (eds) Advances in Geometric Modeling and Processing. GMP 2008. Lecture Notes in Computer Science, vol 4975. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-79246-8_31


f_0 = exactHarmonicOneForm(v, f);

[eta, v_sliced, f_sliced, slice_path] = closedHarmonicOneForm(v, f);

f_0 = [f_0; f_0(slice_path)];

[lambda] = calLambda(v_sliced, f_sliced, f_0, eta);

[energy1, energy2,e2] = calWedgeProductEnergy( v, f, f, f_0, f_sliced, eta);
