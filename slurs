在Linux终端中使用 SLURM，通常是为了提交、监控和管理在计算集群上的作业。以下是一些基本步骤和命令，帮助你在Linux终端中使用SLURM。

### 1. **连接到集群**
首先，你需要通过SSH连接到你的计算集群。如果你已经连接上了，可以跳过这一步。

```bash
ssh your_username@cluster_address
```

### 2. **编写作业脚本**
你需要创建一个作业脚本，这个脚本告诉 SLURM 你的任务需要做什么。通常这是一个 shell 脚本，里面包含了 SLURM 的指令和你的实际工作内容。

```bash
nano my_job_script.sh
```

在脚本中，可以包含如下内容：

```bash
#!/bin/bash
#SBATCH --job-name=my_job          # 任务名称
#SBATCH --output=output.txt        # 标准输出文件
#SBATCH --error=error.txt          # 错误输出文件
#SBATCH --ntasks=1                 # 启动的任务数
#SBATCH --time=01:00:00            # 任务运行的最大时间 (HH:MM:SS)
#SBATCH --partition=normal         # 分区名称 (队列名称)

# 你的实际任务命令
echo "Hello, SLURM!"
```

### 3. **提交作业**
使用 `sbatch` 命令提交你的作业脚本到 SLURM。

```bash
sbatch my_job_script.sh
```

提交后，你会得到一个作业ID，类似于 `Submitted batch job 12345`。

### 4. **查看作业状态**
你可以使用 `squeue` 命令查看当前作业的状态。

```bash
squeue -u your_username
```

这个命令会显示你提交的所有作业的状态。

### 5. **取消作业**
如果你想要取消一个正在运行或者排队的作业，可以使用 `scancel` 命令，并提供作业的ID。

```bash
scancel 12345
```

### 6. **检查作业输出**
作业完成后，SLURM 会将输出保存到你指定的文件中（比如 `output.txt` 和 `error.txt`）。你可以使用 `cat`、`less` 或者其他工具查看这些文件。

```bash
cat output.txt
```

### 7. **常见命令总结**

- **sbatch**：提交作业。
  ```bash
  sbatch my_job_script.sh
  ```
- **squeue**：查看作业队列状态。
  ```bash
  squeue -u your_username
  ```
- **scancel**：取消作业。
  ```bash
  scancel job_id
  ```
- **srun**：直接运行一个命令或者脚本，适合交互式任务。
  ```bash
  srun --ntasks=4 my_program
  ```

这些是 SLURM 在 Linux 终端中的基础操作。根据集群配置和你的需求，你可能还会使用更多的参数和命令。
