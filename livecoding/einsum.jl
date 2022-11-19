#s output_delay = 0.01; prompt_delay=0.2
using OMEinsum
using TropicalNumbers

matmul(A, B) = ein"ij,jk->ik"(A, B);  # define matrix multiplication with einsum

# 定义两个矩阵，第一个参数是类型
A = randn(Float64, 100, 100);
B = randn(Float64, 100, 100);

#s output_delay = 3
@time matmul(A, B);  # matrix multiplication, 1st run

#s output_delay = 0.01
@time matmul(A, B);  # matrix multiplication, 2nd run

# 看看类型推导有没有成功
@code_warntype matmul(A, B)

# 接下来测试 Tropical 代数的矩阵乘法
# 所谓 Tropical 代数，就是把 `*` 操作映射到实数的 `+` 函数
at, bt = Tropical(2.0), Tropical(3.0)
at * bt
# 把 `+` 操作映射到实数的 `max` 函数
at + bt

# 现在把 A 与 B 矩阵的元素类型转成 Tropical 代数类型
At = Tropical.(A);  # `.` means broadcasting, which is similar to Matlab
Bt = Tropical.(B);

#s output_delay = 1.8
@time matmul(At, Bt);  # tropical matrix multiplication, 1st run
#s output_delay = 0.01
@time matmul(At, Bt);  # tropical matrix multiplication, 2nd run

# 看看类型推导有没有成功
@code_warntype matmul(At, Bt)