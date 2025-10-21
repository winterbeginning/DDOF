# trans.dat拟合曲线，用于openfoam13
 通过cantera计算chemkin的trans文件的输运参数，并拟合为在openfoam13中可以使用的形式
## fluidMulticomponentThermo动态库
### ddpolynomial
多项式定义粘度和热导率

### TabulatedFickian

多项式定义质量扩散系数

## ddtransport
读取cantera文件chem.yaml生成质量扩散表格和输运拟合系数

# [原地址](https://github.com/yuchenzh/DDOF)(基于openfoam10的模型)