# Wheat Breeding Optimization with PyMOO

基于 PyMOO 的小麦育种参数优化系统，使用进化算法优化育种程序参数。

## 安装依赖

### Python 依赖

```bash
pip install pymoo rpy2 pyyaml numpy
```

### R 依赖

```r
install.packages(c("AlphaSimR", "yaml"))
```

## evoScript 执行算法流程

### 1. 算法初始化

- 加载配置文件 `config.yml`
- 设置育种参数边界和约束条件
- 初始化 R 环境和 AlphaSimR 库

### 2. 进化算法设置

- 设置使用 NSGA-II 算法

### 3. 目标函数评估

- 对每个个体调用 AlphaSimR 模拟
- 计算遗传增益作为目标函数
- 评估成本约束，超过的给评很差的目标函数值

### 4. 约束检查

- 检查育种参数约束：
  - nCrosses > nDH
  - nPYT > nAYT > nEYT  超过的给评很差的目标函数值

### 5. 进化过程

- 选择优良的个体，执行交叉变异操作，循环执行到到达设置的最大代数

### 6. 结果输出

- 输出每一代进化的参数和遗传增益，最终输出最优的详细参数

## 运行

```bash
conda activate pymoo
python evoScript.py
```
