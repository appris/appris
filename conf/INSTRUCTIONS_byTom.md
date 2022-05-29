# Release 2021_06.v45

Here are guidelines for completing the APPRIS release `2021_06.v45`.

## PART 2

The appristools run is now complete, so it only remains to do quality checks
and, assuming all goes well, upload the new release to the APPRIS server.

## Checking main logs

Check for warnings, errors or exceptions in the main logs at
`/local/appris/tmp/logs/appristools.2021_06.v45.log`,
and the script logs at `/local/appris/tmp/logs/2021_06.v45`.

For Human GENCODE, a small number of warnings are expected relating
to isoform length issues in specific genes (e.g. `MT-ND6`, `MT-CO1`).

## Checking annotation logs

Check the annotations logs as follows.
All commands in this section assume that the
current working directory is `/local/appris/annotations`.

Create a list of annotation log files:
```shell
find 2021_06.v45 -type f -name log > 2021_06.v45_logs.txt
```

### Annotation log file warnings

Create a list of log files containing warning messages:
```shell
while read logfile
do grep -l WARNING $logfile >> 2021_06.v45_warn.txt
done < 2021_06.v45_logs.txt
```

If there are log files containing warnings,
do follow-up searches for warning messages.

Search for `APPRIS::Utils::Logger` warning messages:
```shell
while read logfile
do grep "^WARNING" $logfile >> 2021_06.v45_warn_msgs1.txt
done < 2021_06.v45_warn.txt
```

Warnings may include:
- `WARNING: there are not peptides`: expected in the absence of PROTEO results.
- `WARNING: Only one sequence found.`: most likely due to MEMSAT alignment
  file with a single sequence.
- `WARNING: no relevant Trifid predictions found`:
  expected in the absence of TRIFID results.
- mini-CDS or CDS overshooting the end of a peptide sequence.
- mini-CDS with no residues.

In general, these warnings do not require action.

Search for `APPRIS::Utils::Exception` warning messages:
```shell
while read logfile
do grep -A 1 -- "-- WARNING --" $logfile | grep -v -- "-- WARNING --" >> 2021_06.v45_warn_msgs2.txt
done < 2021_06.v45_warn.txt
```

Warnings may include:
- a residue position exceeding translation length. If the residue position
  exceeds translation length by 1 residue, this warning usually does not
  require action.
- a CDS end overshooting the end of a peptide sequence. If the CDS overshoots
  by 1 residue, this warning usually does not require action.
- a final CDS with no residues. This can happen if the final codon is has 1 or
  2 nucleotides, and in such cases no action is required.

### Annotation log file errors

Create a list of log files containing error messages:
```shell
while read logfile
do grep -il error $logfile >> 2021_06.v45_errs.txt
done < 2021_06.v45_logs.txt
```

If there are log files containing errors,
do a follow-up search for error messages:
```shell
while read logfile
do grep -i error $logfile >> 2021_06.v45_err_msgs.txt
done < 2021_06.v45_errs.txt
```

There is a known error - `ERROR: creating input entity` -
which is believed to be harmless.

### Annotation log file exceptions

Create a list of log files containing exception messages:
```shell
while read logfile
do grep -il exception $logfile >> 2021_06.v45_excs.txt
done < 2021_06.v45_logs.txt
```

If there are log files containing exceptions,
do a follow-up search for exception messages:
```shell
while read logfile
do grep -A 1 -- "-- EXCEPTION --" $logfile | \
       grep -v -- "-- EXCEPTION --" >> 2021_06.v45_exc_msgs.txt
done < 2021_06.v45_excs.txt
```

## Dataset check

Using the script `check_appris_dataset.py`, check
each dataset for consistency as follows.

For Human GENCODE v37:
```shell
python check_appris_dataset.py --double-codon-flag \
    --path /local/appris/data/2021_06.v45/homo_sapiens/e103v45 \
    -o 2021_06.v45_homo_sapiens_e103v45.xlsx
```

For Human GENCODE v38:
```shell
python check_appris_dataset.py --double-codon-flag \
    --path /local/appris/data/2021_06.v45/homo_sapiens/e104v45 \
    -o 2021_06.v45_homo_sapiens_e104v45.xlsx
```

The resulting XLSX files can be inspected for issues with these datasets.
Note that there are currently known issues of moderate severity with
duplicate intervals in Matador3D2 and THUMP.

## Checking breakdown of APPRIS flags

Check the breakdown of Human GENCODE v37 APPRIS flags:
```shell
cd /local/appris/data/2021_06.v45/homo_sapiens/e103v45
awk -F "\t" 'NF == 20 {print $20}' appris_data.appris.txt | sort | uniq -c
```

Check the breakdown of Human GENCODE v38 APPRIS flags:
```shell
cd /local/appris/data/2021_06.v45/homo_sapiens/e104v45
awk -F "\t" 'NF == 20 {print $20}' appris_data.appris.txt | sort | uniq -c
```

## Updating APPRIS server

### Uploading server config file

Upload APPRIS JSON config file from Maliciosa to the APPRIS server:
```shell
rsync -avz /local/appris/conf/config_2021_06.v45.json \
    appris@appris.cnio.es:conf
```

### Update changelog

For updating the changelog for APPRIS release `2021_06.v45`, there is
a Git patch at `/local/appris-aux/patches/appris-2021_06.v45.patch`.

Apply this patch to your working copy of APPRIS, commit the
changes, and push the commit to the APPRIS GitHub repository.
The APPRIS GitHub repository will then be ready for the server update.

### Uploading annotation data

With `/local/appris` having been updated with the latest commits from the
APPRIS GitHub repository (including the updated changelog), upload annotation
data from Maliciosa to the APPRIS server:
```
cd /local/appris
appristools_srv -p 3 -c conf/config_2021_06.v45.json -r 2021_06.v45 \
    -n changelog.md &>tmp/logs/appristools_srv.2021_06.v45.log
```

### Server release

Connect to the APPRIS server:
```shell
ssh appris@appris.cnio.es
```

On the APPRIS server, navigate to the web-server
directory and enter maintenance mode:
```shell
cd /home/appris/server
forever stop $UID && forever start scripts/maintenance.js
```
…where `$UID` is the `forever` process UID of the APPRIS web
server app, which can be obtained by calling `forever list`.

---

At this point, before pulling updates from the APPRIS GitHub repository,
it is typical to stash files that have been modified because they have changes
that are specific to the production server or contain passwords that cannot
be put under source control.

However, on this occasion there are also many files that were modified
on the server to facilitate the transition to HTTPS, so it is necessary
to distinguish these two groups of files and handle them differently.

All changes for the HTTPS update have been incorporated into Git, but
there may be slight differences in some files, so the safest option is
to checkout all these files in order to revert them to the latest commit.

**NB: double-check this command and the list of files before making any changes.**
```
git checkout 2b06b865e208b70ddc18a0c80eddafa29ea0e47e -- \
  INSTALL.md \
  README.md \
  bin/README.md \
  bin/appris_feat_down_ensembl_data \
  code/appris.pl \
  code/src/ensembl/getEComparaAlign.pl \
  code/src/geneset/geneset.pl \
  code/src/spade/spade.pl \
  code/src/trifid/trifid.pl \
  code/src/ucsc/getCDSAlignment.pl \
  code/src/ucsc/getUCSCAlign.pl \
  conf/apprisrc.WS \
  conf/db/apprisdb.sql \
  data/README.md \
  docker/README.md \
  docs/pdoc/main_index.html \
  lib/appris_perllib_v2.1/APPRIS/Analysis/TRIFID.pm \
  lib/appris_perllib_v2.1/APPRIS/Utils/TRIFID.pm \
  scripts/checksum_dataset.pl \
  scripts/insert_appris.pl \
  scripts/run_appris.pl \
  server/README.md \
  server/app/js/app-beta.js \
  server/app/js/app-gold.js \
  server/app/js/app-local.js \
  server/app/js/controllers.js \
  server/app/js/download.js \
  server/app/partials/help/help/methods.html \
  server/app/partials/help/help/services.html \
  server/app/partials/publications.html \
  server/app/templates/navbarTop.tpl.html \
  server/app/templates/tools.tpl.html \
  server/scripts/web-server-beta.js \
  server/scripts/web-server.js \
  ws/.htaccess \
  ws/apidoc/beta.html \
  ws/apidoc/beta/apidoc.json \
  ws/apidoc/beta/exporter.json \
  ws/apidoc/beta/runner.json \
  ws/apidoc/beta/seeker.json \
  ws/apidoc/beta/sequencer.json \
  ws/apidoc/gold.html \
  ws/apidoc/gold/apidoc.json \
  ws/apidoc/gold/exporter.json \
  ws/apidoc/gold/runner.json \
  ws/apidoc/gold/seeker.json \
  ws/apidoc/gold/sequencer.json \
  ws/clients/appris_exporter.pl \
  ws/clients/appris_parser.pl \
  ws/clients/appris_runner.pl \
  ws/functional_isoforms/main.html \
  ws/rest/ebi/muscle_lwp.pl \
  ws/rest/exporter.pl \
  ws/rest/lib/DBRetriever.pm \
  ws/rest/lib/WSRetriever.pm \
  ws/rest/ncbi/web_blast.pl \
  ws/trackHub/docs/APPRIS.html \
  ws/trackHub/docs/CORSAIR.html \
  ws/trackHub/docs/CRASH.html \
  ws/trackHub/docs/FIRESTAR.html \
  ws/trackHub/docs/MATADOR3D.html \
  ws/trackHub/docs/SPADE.html \
  ws/trackHub/docs/THUMP.html \
  ws/trackHub/docs/trackHub.html \
  ws/trackHubPROTEO/docs/PROTEO.html \
  ws/trackHubPROTEO/hub.txt
```

After having checked out the HTTPS-update files,
only the following files should be shown as modified:

- `conf/apprisrc`
- `conf/code/firestar.ini`
- `conf/code/firestar.ini.pro`
- `conf/ws/apprisdb.ini`
- `server/app/gold.html`
- `server/app/index.html`

If so, it should be safe to progress to the next step,
which will involve stashing these specific files.

---

Assuming the working copy of APPRIS is in order,
update the APPRIS server Git repo:
```shell
cd /home/appris
git stash save
git pull origin master
git stash pop
```

Change the current release:
```shell
cd /home/appris/ws/pub
rm current_release && ln -s releases/2021_06.v45 current_release
```

With the JSON config file in the directory `/home/appris/conf`,
make it the new server config file:
```shell
cd /home/appris/conf
rm config.json && ln -s config_2021_06.v45.json config.json
```

Exit maintenance mode:
```shell
cd /home/appris/server
forever stop $UID && forever start scripts/web-server.js
```
…where `$UID` is the `forever` process UID of the APPRIS web
server app, which can be obtained by calling `forever list`.

After exiting maintenance mode, verify that
the APPRIS website is working as expected.

_FIN!_
