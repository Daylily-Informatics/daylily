import os
import sys
import yaml


def t(a):
    os.system(f"touch {a}____x")


t(1)

from snakemake.utils import validate

configfile = "config/config.yaml"
c_yaml_file = open(configfile)
config = yaml.load(c_yaml_file, Loader=yaml.FullLoader)
config["mode"] = "webdav_mount"

import yaml

a_yaml_file = open("config/supporting_files/supporting_dirs.yaml")
parsed_yaml_file = yaml.load(a_yaml_file, Loader=yaml.FullLoader)

config["supporting_dirs"] = parsed_yaml_file
config["logs"] = {"lorm_log": "jemJEMejm.log"}


class SDO(object):
    def __init__(self, config, snakemake, samples):
        self.config = config
        self.snakemake = snakemake
        self.samples = samples
        self.pwd = os.path.abspath(".")
        self.l_root = config["supporting_files"]["link_root"]
        self.l_point = config["supporting_files"]["link_point"]
        self.l_source = config["supporting_files"]["link_source"]

    def _build_local_fs(self):
        for i in sorted(self.config["supporting_dirs"]["root"]):
            print(f"Creating: {i}")
            os.system(
                f"mkdir {i}  > /dev/null 2>&1 || echo '{i} exists' >> {self.config['supporting_dirs']['log']}"
            )

    def local_copy(self):
        self._build_local_fs()
        if True:
            for file_class in self.config["supporting_files"]["files"]:

                for thing_to_copy in self.config["supporting_files"]["files"][
                    file_class
                ]:
                    if self.config["supporting_files"]["files"][file_class][
                        thing_to_copy
                    ]["class"] in ["file_no_suffix_star"]:
                        f_no_suffix_source = (
                            self.config["supporting_files"]["files"][file_class][
                                thing_to_copy
                            ]["name"]
                            .replace("resources/fsx", "/fsx")
                            .split(".")[0]
                        )
                        cmd = f"rsync --times --archive --relative {f_no_suffix_source}* resources/ "
                        print(cmd)
                        # os.system(cmd)
                    elif self.config["supporting_files"]["files"][file_class][
                        thing_to_copy
                    ]["class"] in ["dir_plus"]:
                        d_source = self.config["supporting_files"]["files"][file_class][
                            thing_to_copy
                        ]["name"].replace("resources/fsx", "/fsx")
                        cmd = (
                            f"rsync --times --archive --relative {d_source} resources/ "
                        )
                        print(cmd)
                        # os.system(cmd)
                    elif self.config["supporting_files"]["files"][file_class][
                        thing_to_copy
                    ]["class"] in ["file_only"]:
                        f_source = self.config["supporting_files"]["files"][file_class][
                            thing_to_copy
                        ]["name"].replace("resources/fsx", "/fsx")
                        cmd = (
                            f"rsync --times --archive --relative {f_source} resources/ "
                        )
                        print(cmd)
                        # os.system(cmd)
                    else:
                        raise Exception(
                            f"\n\n\t\tERROR  :: ERROR :: Unexpected local_copy file class: {self.config['supporting_files']['files'][file_class][thing_to_copy]['class']}\n\n................Aborting"
                        )

    # DAYE  Functions
    def local_link(self):
        os.system(f"echo $HOSTNAME__ local_link ")
        os.system("HI!\n\n; sleep 5;")
        os.system(
            f"mkdir {self.l_root} > /dev/null 2>&1 || echo {self.l_root}__exists >> {self.config['logs']['lorm_log']} 2>&1"
        )
        os.system(
            f"sleep 1 "
        )
        #ln -s {self.l_source} {self.l_root} || echo link__{self.l_point}_failed >> {self.config['logs']['lorm_log']} 2>&1"
        

    # DAYE Functions
    def webdav_mount(self):
        os.system(
            'echo """        Warning, you have chosen the webdav mount option for access to supporting files.  This is experimental, and not certain to work if executing on the cloud or cluster."""; sleep 10; '
        )
        c1 = f"mkdir {self.l_root} || echo a >> {self.config['logs']['lorm_log']} 2>&1"
        c2 = f"touch YYYY; mkdir resources/fsxa >> {self.config['logs']['lorm_log']}  2>&1 || echo b >> {self.config['logs']['lorm_log']} 2>&1"
        #c3 = f"etc/lFSXocus.sh {self.l_point} ' -vv ' config/external_tools/.rcc.conf environment/env.sh >> {self.config['logs']['lorm_log']} 2>&1"
        os.system(c1)
        os.system(c2)
        #os.system(c3)


t(2)
# I have no idea how to debug these scripts, can't find the log output or
# start a ipy shell.... so objects wapped in try except it is!
t(3)
print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\t\t\t BEGINNING TO Staging Files\n\n\n")
t(4)

config["supporting_files"]["mode"] = "local_copy"
samples = {}
snakemake = {}
t(5)

sdo = None

try:
    t(6)
    sdo = SDO(snakemake=snakemake, config=config, samples=samples)
    t(7)
except Exception as e:
    print("\n\n\n\n\n\n\t\t\tERROR :: ERROR :: SFO OBJ init failed with", e)
    raise (e)
t(8)
try:
    if "webdav_mount" in [config["supporting_files"]["mode"]]:
        sdo.webdav_mount()
    elif "webdav_copy" in [config["supporting_files"]["mode"]]:
        sdo.wevdav_copy()
    elif "local_link" in [config["supporting_files"]["mode"]]:
        sdo.local_link()
    elif "local_copy" in [config["supporting_files"]["mode"]]:
        sdo.local_copy()
    elif "webdav_direct" in [config["supporting_files"]["mode"]]:
        # Not implemented
        pass
    elif "juicefs_mount" in [config["supporting_files"]["mode"]]:
        # Experimental
        pass
    else:
        print(
            f"\n\n\t\tERROR :: ERROR :: Supporting Data Staging Failed. Unrecognized mode: {snakemake.config['supporting_files']['mode']}"
        )

except Exception as e:
    t(12)
    print(
        "\n\n\n\n\n\n\t\t\tERROR :: ERROR :: Mode Selection Faile for mode:",
        config["supporting_files"]["mode"],
        " ..... \n\n\nand the stack trace",
        e,
    )
    raise (e)
print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\t\t\t DONE Staging Files\n\n\n")
t(13)

#
#  Other options-- dig more deeply into distributed file systems.  I have juicefs working, and it is remarkably performant
#  and fault tolerant.  Like, quite impressive, and designed for the cloud: https://github.com/juicedata/juicefs
#  https://juicefs.com/docs/en/use_juicefs_in_kubernetes.html
