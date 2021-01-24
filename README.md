GSE是基于krylov-decomposition的特征值求解器。krylov分解是Arnoldi分解的进一步推广，它不再需要Arnoldi分解中的如下前提：
1、	U无需保持正交性。
2、	B无需上Hessenberg型。
3、	b不再是em的形式。

GSE的具体算法包含krylov-schur, rational krylov等。

