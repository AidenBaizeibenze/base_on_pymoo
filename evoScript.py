import numpy as np
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2  
from pymoo.optimize import minimize
from pymoo.operators.mutation.pm import PM    
from pymoo.operators.crossover.sbx import SBX 
from pymoo.termination import get_termination
import yaml
import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import tempfile
import random
# import subprocess

class WheatBreedingProblem(Problem):
    def __init__(self, config_path: str = "config.yml"):
        
        self.config = self._load_config(config_path)  
        self.n_factors = self.config["nfactors"]  
        self.name_parameter = self.config["name_parameter"]
        self.is_binary = self.config.get("binary", [False] * self.n_factors)
        self.is_integer = self.config.get("integer", [True] * self.n_factors)
        self.cost_par = self.config.get("cost_par", list(range(1, self.n_factors + 1)))
        self.base_cost_ini = self.config.get("base_cost_ini", [100, 89, 500, 50, 10])
        self.linked_parameter = self.config.get("linked_parameter", False)
        self.which_linked_parameter = self.config.get("which_linked_parameter", [1, 2])
        self.xl, self.xu = self._get_bounds()  
        
        super().__init__(
            n_var=self.n_factors,
            n_obj=1,
            n_ieq_constr=1,
            xl=self.xl,
            xu=self.xu,
        )

        self._setup_r_environment()
    
    def _load_config(self, config_path: str = "config.yml"):
        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
                print(f"config loaded successfully")
                return config  
        except Exception as e:
            print(f"Error config file: {e}")
            return None
    
    def _get_bounds(self):  
        xl = np.array([50, 50, 300, 20, 5, 5])
        xu = np.array([250, 250, 700, 80, 25, 50])
        return xl, xu  
    
    def _check_constraints(self, x):
        ncrosses_dh_constraint = x[0] > x[1]  
        other = x[2] > x[3] and x[3] > x[4]  
        limits = x[2] <= 900 and x[5] <= 50 and x[0] <= 500
        return ncrosses_dh_constraint and other and limits

    def _evaluate(self, X, out, *args, **kwargs):
        n_individuals = X.shape[0]
        
        F = np.zeros((n_individuals, self.n_obj))
        G = np.zeros((n_individuals, self.n_ieq_constr))
        
        for i in range(n_individuals):
            x = X[i, :]
            if not self._check_constraints(x):
                F[i, 0] = 1000  
                print(f"个体{i}违反: {x}")
            else:
                F[i, 0] = self.alphaSimR(x)  
            
            G[i, 0] = self._cost_fun(x)  
        
        out["F"] = F
        out["G"] = G

    #测试用
    # def _objective_fun(self, x: np.ndarray) -> float:  
    #     # 替换为AlphaSimR
    #     genetic_gain = (
    #         x[0] * 0.1 + x[1] * 0.05 + x[2] * 0.08 + 
    #         x[3] * 0.12 + x[4] * 0.15 
    #     )
    #     return -genetic_gain
    
    def _setup_r_environment(self):
        try:
            robjects.r('library(AlphaSimR)')
            print("AlphaSimR_loaded")
        except Exception as e:
            print(f"R_failed: {e}")
    
    def alphaSimR(self, x: np.ndarray) -> float:  
        try:
            
            return self._call_alphasimr(x)
        except Exception as e:
            print(f"AlphaSimR_failed: {e}")
            
    
    def _call_alphasimr(self, x: np.ndarray) -> float:
        try:
            with tempfile.NamedTemporaryFile(suffix='.RData', delete=False) as tmp_file:
                file = tmp_file.name
            
            rep = 1
            random_seed = random.randint(10000, 99999)  
            nCrosses = int(x[0])
            nDH = int(x[1])
            nPYT = int(x[2])
            nAYT = int(x[3])
            nEYT = int(x[4])
            newParents_replace = int(x[5])
            
            #call R script
            r_script = f"""
            
            args <- c('{file}', {rep}, {random_seed}, {nCrosses}, {nDH}, {nPYT}, {nAYT}, {nEYT}, {newParents_replace})
            
            outcome <- args[1]
            rep <- as.integer(args[2])
            randomSeed <- as.integer(args[3]) * rep

            nCrosses <- as.numeric(args[4])
            nDH <- as.numeric(args[5])
            nPYT <- as.numeric(args[6])
            nAYT <- as.numeric(args[7])
            nEYT <- as.numeric(args[8])
            newParents_replace <- as.numeric(args[9])
            
            print(paste("randomSeed:", randomSeed))
            
            results <- NULL
            library(AlphaSimR)
            
            assign("commandArgs", function(trailingOnly=FALSE) {{
                return(args)
            }}, envir = .GlobalEnv)
            
            source('simuScript.r')
            save(list=c("results"), file = outcome)
            """
            
            robjects.r(r_script)
            robjects.r(f"load('{file}')")
            results = robjects.r('results')
            genetic_gain = float(robjects.r(f"results[1,9]")[0])
            print(f"遗传增益: {genetic_gain}")
            os.unlink(file)
            
            return -genetic_gain
            
        except Exception as e:
            print(f"AlphaSimR调用错误: {e}")
            return 0


    def _cost_fun(self, x):
        current_cost = self._caculate_cost(x)
        base_cost = self._caculate_cost(self.base_cost_ini)
        cost_difference = abs(current_cost - base_cost)  
        cost_constraint = cost_difference - 5000
        return cost_constraint
    
    def _caculate_cost(self, x):
        cost = (
            x[0] * 30 +  # nCrosses 
            x[0] * 30 +  # F1 
            x[1] * x[0] * 30 +  # DH 
            x[1] * x[0] * 15 +  # GS
            x[2] * 5 * 20 +  # PYT 
            x[3] * 15 * 50 +  # AYT 
            x[4] * 20 * 50    # EYT 
        )
        return cost
    
    def _test(self):
        print("测试定义")
        test_x = np.array([100, 80, 500, 50, 10, 20])
        
        objective = self.alphaSimR(test_x)  
        print(f"目标函数值: {objective}")
        # constraint = self._cost_fun(test_x)  
        # print(f"约束函数值: {constraint}")
        # cost = self._caculate_cost(test_x)  
        # print(f"成本: {cost}")
        
        print("done！")


class WheatBreedingAlgorithm:
    def __init__(self, config_path: str = "config.yml"):
        self.problem = WheatBreedingProblem(config_path)
        self.mut_parent = self.problem.config["mut_parent"]
        self.offspring = self.problem.config["mut_offspring"]
        self.results = None

    def optimize(self, population_size: int = 1000, n_generations: int = 150):  
        algorithm = NSGA2(
            pop_size=population_size,
            crossover=SBX(prob=0.9, eta=15),
            mutation=PM(prob=self.offspring, eta=20),
            eliminate_duplicates=True  
        )

        termination = get_termination("n_gen", n_generations)  
        res = minimize(
            self.problem,
            algorithm,
            termination,
            seed=1,
            save_history=True,
            verbose=True,
        )
        self.results = res
        return res
    
    def print_results(self):
        if self.results is None:
            print("没有优化结果")
            return
        
        if len(self.results.F) == 1:
            best_solution = self.results.X  
            best_objective = self.results.F[0]
        else:
            best_idx = np.argmin(self.results.F)
            best_solution = self.results.X[best_idx]
            best_objective = self.results.F[best_idx]
        
        print(f"\n" + "="*5 + "优化结果" + "="*5)
        print(f"最优解: {best_solution}")
        print(f"最优目标值: {best_objective}")
        
        print(f"\n" + "="*5 + "参数结果" + "="*5)
        for i, name in enumerate(self.problem.name_parameter):
            if i < len(best_solution):
                print(f"{name}: {int(best_solution[i])}")





if __name__ == "__main__":
    problem = WheatBreedingProblem()
    problem._test()  
    
    print("\n" + "="*50)

    optimizer = WheatBreedingAlgorithm()
    results = optimizer.optimize(population_size=10, n_generations=10)  
    optimizer.print_results()