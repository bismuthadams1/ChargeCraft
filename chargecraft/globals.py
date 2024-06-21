import os  
import psutil

class GlobalConfig:
    num_threads_per_core: int = 2  # Example default value
    memory_allocation: float = 2e+9  # Example default value
    
    @staticmethod
    def total_threads():
        num_cores = int(os.getenv('SLURM_CPUS_PER_TASK', 8))  # Corrected typo here
        return num_cores * GlobalConfig.num_threads_per_core
    
    @staticmethod
    def memory():
        return GlobalConfig.memory_allocation
    

def log_memory_usage():
    process = psutil.Process()
    print(f"Memory Usage: {process.memory_info().rss / 1024 ** 2} MB", flush=False)