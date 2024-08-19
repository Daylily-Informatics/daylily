# Cost Tags 
- The default project cost tags are `daylily-dev`, and `daylily-prd` may be set as well.
- To change the tag, change the `DAY_PROJECT` environment variable after running `dyiniy` and prior to running `dy-r` :
```bash
export DAY_PROJECT="daylily-dev"
```
- Only allowed project tag values will be allowed to run via slurm (tags have no effect on locally run jobs). The allowed tag list should be edited on the headnode in the following file: ` /opt/slurm/etc/projects_list.conf`. Add new project tags to the daylily user to enable them.
