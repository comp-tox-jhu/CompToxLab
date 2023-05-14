## Unix Cheatsheet

## File Commands

| Command | Description |
|---|---|
| ls |  directory listing| 
| ls -al |  formatted listing with hidden files| 
| cd dir |  change directory to dir| 
| cd |  change to home| 
| pwd | show current directory| 
| mkdir dir |  create a directory dir| 
| rm file |  delete file| 
| rm -r dir |  delete directory dir| 
| rm -f file |  force remove file| 
| rm -rf dir |  force remove directory dir *| 
| cp file1 file2 |  copy file1 to file2| 
| cp -r dir1 dir2 |  copy dir1 to dir2; create dir2 if it doesn't exist| 
| mv file1 file2 |  rename or move file1 to file2. if file2 is an existing directory, moves file1 into directory file2| 
| ln -s file link |  create symbolic link link to file| 
| touch file |  create or update file| 
| cat > file |  places standard input into file| 
| more file |  output the contents of file| 
| head file |  output the first 10 lines of file| 
| tail file |  output the last 10 lines of file| 
| tail -f file|  output the contents of file as it grows, starting with the last 10 lines Process Management| 

## Process Management

| Command | Description |
|---|---|
| ps |  display your currently active processes| 
| top |  display all running processes| 
| kill pid |  kill process id pid| 
| killall proc |  kill all processes named proc *| 
| bg |  lists stopped or background jobs; resume a stopped job in the background| 
| fg |  brings the most recent job to foreground| 
| fg n |  brings job n to the foreground| 


## File Permissions

chmod guo+ file – change the permissions of file/folder for the group (g), user (u), other (o)

- read (r)
- write (w)
- execute (x)

Examples:

- chmod gu+rwx file – read, write, execute for everyone in the group and user (file)
- chmod -R gu+rwx folder – read, write, execute for everyone in the group and user (folder)
- chmod g-w - remove write for everyone in the group


## SSH

| Command | Description |
|---|---|
| ssh user@host |  connect to host as user| 
| ssh -p port user@host |  connect to host on port, port as user| 
| ssh-copy-id user@host |  add your key to host for user to enable a keyed or passwordless login| 


## Searching

| Command | Description |
|---|---|
| grep pattern files |  search for pattern in files| 
| grep -r pattern dir |  search recursively for pattern in dir| 
| command | grep pattern |  search for pattern in theoutput of command| 
| locate file |  find all instances of file| 


## System Information

| Command | Description |
|---|---|
| date |  | show the current date and time| 
| cal |  show this month's calendar| 
| uptime |  show current uptime| 
| w |  display who is online| 
| whoami |  who you are logged in as| 
finger user | – display information about user| 
| uname -a | show kernel information| 
| cat /proc/cpuinfo |  cpu information| 
| cat /proc/meminfo |  memory information| 
| man command |  show the manual for command| 
| df |  show disk usage| 
| du |  show directory space usage| 
| free |  show memory and swap usage| 
| whereis app |  show possible locations of app| 
| which app |  show which app will be run by default| 


## Compression

| Command | Description |
|---|---|
| tar cf file.tar files |  create a tar named file.tar containing files| 
| tar xf file.tar |  extract the files from file.tar| 
| tar czf file.tar.gz files |  create a tar with Gzip compression| 
| tar xzf file.tar.gz |  extract a tar using Gzip| 
| tar cjf file.tar.bz2 |  create a tar with Bzip2 compression| 
| tar xjf file.tar.bz2 |  extract a tar using Bzip2| 
| gzip file |  compresses file and renames it to file.gz| 
| gzip -d file.gz |  decompresses file.gz back to file| 


## Network 

| Command | Description |
|---|---|
| ping host |  ping host and output results| 
| whois domain |  get whois information for domain| 
| dig domain |  get DNS information for domain| 
| dig -x host |  reverse lookup host| 
| wget file |  download file| 
| wget -c file |  continue a stopped download| 

## Shortcuts

| Command | Description |
|---|---|
| Ctrl+C |  halts the current command| 
| Ctrl+Z |  stops the current command, resume with fg in the foreground or bg in the background| 
| Ctrl+D |  log out of current session, similar to exit| 
| Ctrl+W |  erases one word in the current line| 
| Ctrl+U |  erases the whole line| 
| Ctrl+R |  type to bring up a recent command| 
| !! |  repeats the last command| 
| exit |  log out of current session| 


## References

1. [Unix/Linux Command Reference](https://upload.wikimedia.org/wikipedia/commons/7/79/Unix_command_cheatsheet.pdf)
