新网格模拟：Jell、载荷核心谱、因子协方差与方向恢复审计
数据只读；无重新拟合。全部10796个compatible有效估计均纳入，4个缺失另由完整性审计处理。

最主要结论
1. 外部J没有显示趋零；按实际分解basis匹配后，J谱向各basis的正值总体谱靠近。
   原五列basis J_min均值T100/200/400/800为 .017206/.015836/.015270/.014994，
   对应总体 .014665831。最后一个T估计/参考均值比1.0224。
   其余实际出现basis在T800：
   1_2_3_5_6：n7，J_min=.028221，参考=.027171；
   1_2_4_5_6：n23，J_min=.038080，参考=.037092；
   1_3_4_5_6：n16，J_min=.018223，参考=.016993。
   原basis n854，替代basis合计46/900。总体未出现的合法basis不从样本推断。
2. 外部混合所有basis的J_min均值为 .019371/.017479/.016406/.015744；
   对应实际basis人口参考均值为 .016289/.016101/.015695/.015378。
   该均值下降同时反映估计误差减小及basis组成改变，不能解释为趋于0。
   全J谱orig平均绝对误差 .01543→.003518，nor .01306→.003092。
3. 边界J_min均值 .019290/.014366/.010785/.007977，向0下降；
   内部 .088442/.087701/.087360/.087186，保持正值。
   nor对应内部 .062734/.062147/.061893/.061764；
   边界 .014087/.010489/.007869/.005818；
   外部 .026312/.024043/.022666/.021762。
4. T800各固定(p,q)的三类Jnor最小值范围均不重叠。
   但把不同(p,q)混合后，边界最大.014069略超过外部最小.013892，不能称合并范围完全分离。

Ar、Br、Z、C（Omega）
5. 估计Ar/Br在1e-10数值阈值下全为满列秩：内部8，另两组5。
   这不代表总体rank变化：真实Ar/Br rank分别内部4、另两组3。
   底部小奇异值随T减小；内部Ar底4个范数 .02249→.007697，
   边界Ar底2个 .01617→.005607，外部Ar底2个 .01908→.005883。
6. C（Omega的因子矩阵）总体内部rank26/28，边界/外部9/10。
   三种C都有趋零尾谱：内部底2个范数 .023546→.008052；
   边界底1个 .014500→.005421；外部底1个 .006822→.001828。
   因此Omega欠秩本身不是内外区分标准；内部仍有正J谱。
7. Z始终满列秩。原basis总体最小奇异值内部.406020，边界.217039，外部.151646。
   外部原basis条件更弱，可解释其有限样本几何残差较难恢复。

Mx与信号
8. 全部估计Mx特征值保持正值。全样本最小Mx特征值约 .02521（外部），
   不能把替代basis中的小Mx特征值直接认作弱因子：Mx随可逆basis变换而改变。
   物理信号矩阵F=ZL的最小奇异值/sqrt(pq)，外部全样本最低T100 .14172、T800 .15568。
   没有旧坏例子Mx~1e-10那样的近丢失因子退化。
   T800均值：内部.55273，边界.30432，外部.21274。
9. 命名必须区分：J/Ar/Br/Z/C表nor = s_i/s_max；
   fit_covariance_spectra表signal_singular_value_nor = s_i(F)/sqrt(pq)，不是除以最大奇异值。
   Mx_eigenvalue本来就是对称Mx特征值，不平方。
   J_gram_eigenvalues另表 = J_orig^2/J_nor^2；已逐行验证平方关系误差<1e-12。

方向恢复
10. 每个估计列对已知真方向进行不重复匹配，外部允许6选5，不按是否回原basis删除。
    外部每列匹配余弦均>=.98的估计数为779/900、860/900、886/900、895/900。
    边界891/900后900/900；内部T100为896/896，其后900/900。
    外部T800最差列余弦.94186，集中小p/q替代basis；不应称全部900都准确恢复。
11. 外部额外物理方向到估计载荷空间的投影残差均值 .09721→.03344，
    投影E1均值 .08078→.02886。这些是使用模拟真值的oracle诊断，不是实际数据检验。

p、q效应及旧实验比较
12. T800从(20,20)到(80,80)的空间误差：
    内部 .01060→.00537；边界 .01770→.00892；外部 .02385→.01210。
    (20,80)与(80,20)的均值几乎对称。具体9组见metric_spectra_summary_by_pq_T.csv。
13. 新旧40x40内部空间误差T200/400/800几乎完全相同，
    但新strict收敛标志0%/3%/16%，旧91%/100%/100%，说明收敛flag低不等于参数恢复差。
    边界T100空间误差从旧.10839降到新.03641；旧小T差解明显减少。
    外部新旧分支组成不同：T800新extra选入2%，旧31%，导致J平均值不能不分basis直接对比。

主要输出
all_J_singular_values_mean_variance.csv：每个p/q/T/regime/ell/奇异值序号的orig/nor及Gram，均值和样本方差。
J_spectrum_by_actual_basis.csv：按实际basis与对应真方向ID拆分的全部谱。
Jell_minimum_singular_values.csv：每个ell最小奇异值的均值/方差及参考。
exterior_J_reference_by_seed.csv / exterior_J_reference_summary.csv：逐估计及各p/q的参考比较。
exterior_J_by_actual_basis_summary.csv：实际basis汇总。
Ar_Br_Z_C_all_spectra_summary.csv：四类矩阵全部谱。
matrix_rank_and_tail_summary.csv / matrix_rank_and_tail_pooled.csv：数值秩及真秩以下的尾谱。
Mx_eigenvalues_and_signal_svd_summary.csv：全部Mx特征值、信号奇异值。
direction_recovery_by_pq_T.csv / exterior_basis_counts.csv：方向匹配、原/替代basis计数。
extra_true_direction_projection_summary.csv：额外真方向投影诊断。
metric_spectra_summary_by_pq_T.csv：所有组合的主谱指标、误差与匹配mean/var。
old_vs_new_40x40_spectra_summary.csv：旧/新同维度比较。

所有方差均采用分母n-1。不建议把不同basis下的同ell自动当作同一个物理方向；
比较物理分支时使用assigned_true_direction_id和population_basis_key。
