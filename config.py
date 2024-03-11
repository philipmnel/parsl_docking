from parsl.config import Config
from parsl.executors.threads import ThreadPoolExecutor
from parsl.utils import get_all_checkpoints
from parsl.executors import HighThroughputExecutor
from parsl.providers import SlurmProvider
from parsl.launchers import SrunLauncher

tconfig = Config(executors=[ThreadPoolExecutor(
    max_threads=1, label='local_threads'
)],
#checkpoint_mode='task_exit')
#checkpoint_files = get_all_checkpoints())
)

ht_config = Config(
	executors=[
		HighThroughputExecutor(
			label="htex",
			worker_debug=False,
			cores_per_worker=8.0,
			max_workers=1,
			provider=SlurmProvider(
				partition="hive",
				account="hive-cs207",
				nodes_per_block=1,
				mem_per_node=32,
				init_blocks=1,
				max_blocks=16,
				scheduler_options="",
				cmd_timeout=60,
				walltime="01:00:00",
				launcher=SrunLauncher(),
				worker_init="module load anaconda3; conda activate vina",
			),
		)
	],
)
