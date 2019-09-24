# Note: this program must be run using bash.  It will not run in fish.

function BoxMount() {

    #local LOGIN_FILE='/home/rweaton/Documents/.secret'
    local LOGIN_FILE='/home/thugwithyoyo/Documents/.secret'

    login=($(< $LOGIN_FILE))
    local USER=${login[0]}
    local PASSWD=${login[1]}

    curlftpfs ftp.box.com ~/BoxFTP -o user=$USER:$PASSWD,allow_other

    #rsync -r ~/Dropbox/CNPRC/ "/home/rweaton/BoxFTP/MoxonLab/Projects/NHP Calcium Imaging/SharedDropboxSync/CNPRC"

    #sudo fusermount -u ~/BoxFTP

}

function GeneratePathsList() {
	
    local ScriptDir=$(pwd)
    local TargetDir="/home/thugwithyoyo/BoxFTP/MoxonLab/Projects/NHP Calcium Imaging/NHP_CalciumImaging/PlexonData/RetrackingData"
    cd "${TargetDir}"
    find . -name '*.csv' > ~/Desktop/PathsToCSVs.txt
    cd "${ScriptDir}"

}

BoxMount
#GeneratePathsList

#TargetDir="~/BoxFTP/MoxonLab/Projects/NHP Calcium Imaging/NHP_CalciumImaging/PlexonData/RetrackingData"

#umount -f ~/BoxFTP
