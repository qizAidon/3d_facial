 %1: 这是一份说明文档
 %2: Author: Zheng Qi
 %3: Date: 11/17/2016
 
 %4: 该文件夹包含./files子文件夹和./spin_imgs子文件夹，以及若干.m文件，用于实现
%论文中所描述的算法。

 %5: ./file文件夹中包含所需要的./pts文件和./abs文件，两种文件包含的具体信息见
%meshSpin.m的函数说明。

 %6：./spin_imgs包含了实验中用到的所有50个样本的spin image特征，因为每一个mesh
%样本上存在两万~三万个顶点，计算该特征耗时很长。实际应用中也可以只计算所需顶点的
%spin image特征，这可以作为后期的一个修改。

 %7：main.m中按照算法流程标注了大致步骤，是非常粗略的版本，因为在做PCA时遇到了
%问题还没有解决（论文中有关于该问题的详细描述），运行时会出现
%“警告: 矩阵接近奇异值，或者缩放错误...”
%因此，该main.m还没有进行模块化整体以及局部优化等。

 %8：calcSpinImages.m用于计算输入mesh上每个顶点的spin image特征，是作者或者第三方
%对该特征的MATLAB实现，我根据需要做了一些小修改，也补充了函数说明。

 %9：extract_abs.m用于读取.abs文件中的信息，包括flag和三维坐标。
 
 %10：read_shape.m用于读取.pts文件中关键点的二维坐标。
 
 %11：meshSpin.m利用前面函数获取每一个输入mesh的spin image特征，三维坐标，关键点
%的三维坐标以及关键点的spin image.注意，对于测试样本，其.pts输入参数应为空，因为
%测试样本关键点的位置是未知的，这一点在meshSpin.m的函数说明里面有。而对于训练样本,
%应输入.pts文件和对应的.abs文件。