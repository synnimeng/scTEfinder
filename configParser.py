import pandas as pd
from types import SimpleNamespace
import os


Option_Platform = ('10x3-v1','10x3-v1.r3',
                   '10x3-v2','10x3-v3',
                   '10x5','10x5-v1','10x5-v2')

# 判断单元格为空
def isNullCell(v):
    return pd.isna(v) or (str.strip(str(v)) == "")


# region |- 工具配置项 -|
class Tool():
    def __init__(self, path) -> None:
        table = pd.read_csv(path, sep='\t', index_col=0)
        self.STAR = table.loc['STAR', 'path']
        self.seqkit = table.loc['seqkit', 'path']
        self.Rscript = table.loc['Rscript', 'path']
        self.python = table.loc['python', 'path']
        self.bcsubset = table.loc['bcsubset', 'path']
        self.scTE = table.loc['scTE', 'path']



class Job():
    def __init__(self, info:pd.DataFrame):
        # "DataDir" "Postfix" "Platform" "RefGenome" "OutDir" "Memory" "Thread" "Description"
        self.Job = info['Job']
        self.DataDir = info['DataDir']
        self.Prefix = info.get('Prefix', '')
        self.InPrefix = os.path.join(self.DataDir, self.Prefix)
        self.Postfix = info['Postfix']
        self.Platform = info['Platform'].lower()
        assert self.Platform in Option_Platform, f"Platform {self.Platform} Not In Option: {Option_Platform}"
        self.Tissue = info['Tissue'].lower().replace(',', ' ')  # , 转 空白, 符合nargs 输入
        self.R1FILE = f"{self.InPrefix}_R1.{self.Postfix}" if isNullCell(info.get("R1", "")) \
            else f'{self.DataDir}/{info["R1"]}'
        self.R2FILE = f"{self.InPrefix}_R2.{self.Postfix}" if isNullCell(info.get("R2", "")) \
            else f'{self.DataDir}/{info["R2"]}'
        # 仅部分10x3-v1平台可用
        self.RUFILE = None
        if self.Platform.endswith('r3'):
            self.RUFILE = f"{self.InPrefix}_RU.{self.Postfix}" if isNullCell(info.get("RU", "")) \
                else f'{self.DataDir}/{info["RU"]}'
        self.RefGenome = info['RefGenome']
        self.RefGenomeIndex = info['RefGenomeIndex']
        self.CellTypistModel = info['CellTypistModel']
        self.OutDir = info['OutDir']
        self.Memory = int(info['Memory'])
        self.Thread = int(info['Thread'])
        self.CondaEnv = info['CondaEnv']
        self.Description = info['Description']



class JobQueue():
    def __init__(self, path, path_default=None):
        self.workDir = os.getcwd()
        self.defaultSets = pd.read_csv(path_default, sep='\t').fillna("") if path_default else None
        self.schedule = pd.read_csv(path, sep='\t').fillna("")
        assert self.schedule["Job"].duplicated().sum()==0, "Job Name exist Duplicate Value!\nPlease check and rename it."
        # 补充默认值
        for k in ("DataDir", "Postfix", "Tissue", "Platform", "RefGenome", "RefGenomeIndex", "CellTypistModel", 
                  "OutDir", "Memory", "Thread", "CondaEnv", "Description"):
            if k not in self.schedule.keys(): self.schedule[k] = self.defaultSets[k].values[0]
            # 填充空值
            self.schedule[k] = self.schedule[k].apply(lambda x: self.defaultSets[k].values[0] if isNullCell(x) else x)
        # 修改相对路径 > 绝对路径
        for k in ("DataDir", "OutDir", "RefGenome", "RefGenomeIndex"):
            self.schedule[k] = self.schedule[k].apply(lambda x: os.path.join(self.workDir, x))
        self.jobs = {
            info["Job"]: Job(info) for i, info in self.schedule.iterrows()
        }
    def __getitem__(self, key):
        return self.schedule[key].tolist()
    def getjob(self, job, safe=True):
        if safe:
            assert job in self.jobs, f"Job {job} Not Exist! When getjob()"
        return self.jobs.get(job, None)    # 默认值

