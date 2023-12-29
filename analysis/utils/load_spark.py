import os
import random
from itertools import count
import hail as hl
from uuid import uuid4
import platform
from pyspark.sql import SparkSession


WD = None
SC = None
RESOURCES_PATH = None

if 'aws' in platform.release():
    builder = (
        SparkSession
        .builder
        .config('spark.driver.memory', '60G')
        .config('spark.executor.memory', '60G')
        .config('spark.ui.showConsoleProgress', 'false')
        .enableHiveSupport()
    )
    SC = builder.getOrCreate()

    hl_init_kwargs = {
        'sc': SC.sparkContext,
        'default_reference': 'GRCh38'
    }
elif '.cyf.' in platform.release():
    master = None
    spark_conf = {
        'spark.driver.memory': '40G',
        'spark.executor.memory': '80G',
        'spark.rpc.message.maxSize': '256',
    }

    localfs_path = os.environ.get('SCRATCH_LOCAL') + '/'
    scratch_path = os.environ.get('SCRATCH') + '/'
    hail_log_uuid = str(uuid4())

    # ares
    if platform.node().startswith('ac'):
        WD = '/net/pr2/projects/plgrid/plggneuromol/ifpan-gosborcz-ukb/'
        RESOURCES_PATH = '/net/pr2/projects/plgrid/plggneuromol/resources/'
        os.environ['_JAVA_OPTIONS'] = f'-Djava.io.tmpdir={localfs_path}'
    # prometheus
    elif platform.node().startswith('p'):
        WD = '/net/archive/groups/plggneuromol/ifpan-gosborcz-ukb/'
        RESOURCES_PATH = '/net/archive/groups/plggneuromol/resources/'
        spark_master_host = os.environ.get('SPARK_MASTER_HOST')
        spark_master_port = os.environ.get('SPARK_MASTER_PORT')
        master = f'spark://{spark_master_host}:{spark_master_port}'

    hl_init_kwargs = {
        'master': master,
        'tmp_dir': os.path.join(scratch_path, 'hail-tmpdir'),
        'default_reference': 'GRCh38',
        'spark_conf': spark_conf,
        'log': os.path.join(scratch_path, f'slurm-log/hail-{hail_log_uuid}.log'),
        'local_tmpdir': os.path.join(localfs_path, 'hail-local-tmpdir')
    }
else:
    raise ValueError(f'Platform release {platform.release()} unknown.')


def hl_init(**kwargs):
    hl_init_kwargs.update(kwargs)
    hl.init(**hl_init_kwargs)


def tmpdir_path_iter(prefix=None):
    if int(os.getenv('SLURM_NNODES')) > 1 or os.getenv('HAIL_CHECKPOINT_ENV'):
        tmp_path = os.path.join(scratch_path, 'tmp/')
    else:
        tmp_path = os.path.join(localfs_path, 'tmp/')
    os.makedirs(tmp_path, exist_ok=True)
    if prefix is None:
        prefix = f"{random.randrange(16 ** 4):04x}"
    counter = count()
    while True:
        path = os.path.join(
            tmp_path,
            f"{prefix}-{os.getenv('SLURM_JOBID')}-{str(next(counter))}"
        )
        if os.getenv('HAIL_CHECKPOINT_ENV'):
            print(f'HAIL_CHECKPOINT_PATH: {path}', flush=True)
        yield path

wd = WD
