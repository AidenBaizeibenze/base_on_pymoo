# Wheat Breeding Optimization with PyMOO

基于 PyMOO 的小麦育种参数优化系统，使用进化算法优化育种程序参数。

## 功能特点

- 使用 NSGA-II 进化算法优化育种参数
- 集成 AlphaSimR 进行遗传模拟
- 支持约束优化和成本控制
- 可配置的育种参数和约束条件

## 文件结构

```
base_on_pymoo/
├── evoScript.py      # 主程序文件
├── simuScript.r      # R 模拟脚本
├── config.yml        # 配置文件
├── environment.yml   # 环境依赖
└── README.md         # 说明文档
```

## 安装依赖

### Python 依赖
```bash
pip install pymoo rpy2 pyyaml numpy
```

### R 依赖
```r
install.packages(c("AlphaSimR", "yaml"))
```

## 使用方法

1. 配置 `config.yml` 文件
2. 运行优化程序：
```python
python evoScript.py
```

## 参数说明

- `nCrosses`: 杂交数量
- `nDH`: 双单倍体数量  
- `nPYT`: 初步产量试验数量
- `nAYT`: 高级产量试验数量
- `nEYT`: 精英产量试验数量
- `newParents_replace`: 新亲本替换数量

## 约束条件

- nCrosses > nDH
- nPYT > nAYT > nEYT
- 成本约束：与基准成本差异不超过 5000

## 许可证

MIT License
