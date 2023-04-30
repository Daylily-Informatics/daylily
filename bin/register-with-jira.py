#!/usr/bin/env python3
"""Script to take details from a snakemake recipe and register with Jira once complete"""

import os
import sys

jira = str(sys.argv[1]).replace('"', "-")
details = (
    os.popen("echo $DAY_PROFILE $DAY_BIOME").readline().rstrip()
)  # str(sys.argv[2]).replace(' ', "-").replace('"','')
EX = str(sys.argv[3]).replace('"', "-")
RU = str(sys.argv[4]).replace('"', "-")
anapipe = str(sys.argv[5]).replace('"', "-")
jobstate = str(sys.argv[6]).replace('"', "-")
url = str(sys.argv[7]).replace('"', "-")
nsamp = str(sys.argv[8]).replace('"', "-")
c1 = str(sys.argv[9]).replace('"', "-")
c2 = str(sys.argv[10]).replace('"', "-")
c3 = str(sys.argv[11]).replace('"', "-")

emojiia = " (flagoff) (flagoff) (flagoff) "
if jobstate == "success":
    emojiia = " (/) "
elif jobstate in ["fail", "failed"]:
    emojiia = " (x) "

datadirs = os.popen("cwd").readline().rstrip()

cmd = (
    '''curl -X POST -H 'Content-type: application/json' --data '{"issues":["'''
    + jira
    + '''"],"details":"'''
    + details
    + '''","EX":"'''
    + EX
    + '''","Run":"'''
    + RU
    + '''","anapipe":"'''
    + anapipe
    + '''","jobstate":"'''
    + jobstate
    + '''","emojiia":"'''
    + emojiia
    + '''","datadirs":"'''
    + datadirs
    + '''","nsamp":"'''
    + nsamp
    + '''","c1":"'''
    + c1
    + '''","c2":"'''
    + c2
    + '''","c3":"'''
    + c3
    + '''","url":"'''
    + url
    + """"}' "https://automation.atlassian.com/pro/hooks/c05058e99fb6992714547689d81d5e310c1cc797" """
)

# cmd = """curl -X POST -H 'Content-type: application/json' --data '{7}"issues":["TEST-1"], "data": {7}"epic_link":"SUP-48", "sample": "{0}", "analysis_type" : "{1}{2}", "results_url" : "{3}", "time_to_execute": "{4}", "other_stuff":"{5}", "notify_users": "{6}" {8} {8}' "http"""".format(    sample,    analysis_type,    treatment,    url,    ex_time,    other_stuff,  jira_users_to_tag,    "{",    "}",)
print(cmd)
res = os.popen(cmd).readlines()

print(res)
