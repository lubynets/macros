alias gsi='ssh -X lubynets@lxi001.gsi.de'
alias mntgsi='sshfs lubynets@lxi001.gsi.de:/u/lubynets/ /home/oleksii/Mount/gsi/'
alias mntlustre2steps='sshfs lubynets@lxi001.gsi.de:/u/lubynets/Mount/Lustre/ /home/oleksii/Mount/Lustre/'
alias mntlustre="sshfs -o ssh_command='ssh -J lubynets@lxi001.gsi.de' lubynets@lustre.hpc.gsi.de:/lustre/cbm/users/lubynets/ /home/oleksii/Mount/Lustre/"
# alias mntlustre="sshfs -o ssh_command='ssh -J lubynets@lxi001.gsi.de' lubynets@vae22.hpc.gsi.de:/lustre/cbm/users/lubynets/ /home/oleksii/Mount/Lustre/"
# alias mntlustre="sshfs -o ssh_command='ssh -J lubynets@lxi001.gsi.de' lubynets@virgo-centos7.hpc.gsi.de:/lustre/cbm/users/lubynets/ /home/oleksii/Mount/Lustre/"
alias mmg="mntgsi && mntlustre && gsi"
alias umntgsi='sudo umount ~/Mount/gsi'
alias umntlustre='sudo umount ~/Mount/Lustre'
alias kssh='killall -9 ssh'

alias cern='ssh -X ollubyne@lxplus.cern.ch'
alias mntcern='sshfs ollubyne@lxplus.cern.ch:/afs/cern.ch/user/o/ollubyne/ /home/user/Mount/cern/'

alias mji='make -j3 install'
alias mj='make -j3'
alias mi='make install'
alias m='make'
alias cd1='cd ..'
alias cd2='cd ../..'
alias cd3='cd ../../..'
alias cd4='cd ../../../..'
alias cd5='cd ../../../../..'
alias cd6='cd ../../../../../..'
alias cd7='cd ../../../../../../..'

alias gst='git status'
alias gpl='git pull'
alias glo='git log'
alias gdf='git diff'
alias gch='git checkout'
alias grv='git remote -v'
alias kt='killall -9 telegram'
alias clang-format='~/.clang-format'
alias pdf2png='pdftoppm -png -cropbox -singlefile'
alias pdf2pngall='pdftoppm -png -cropbox'
alias cpdf='~/.cpdf' # removes pictures from pdf: cpdf -draft input.pdf -o output.pdf

alias AT='source /home/oleksii/soft/AnalysisTree/install_master/bin/AnalysisTreeConfig.sh'
alias qntoolsconfig='source /home/oleksii/soft/.qntoolsconfig.sh'
alias flowdrawconfig='source /home/oleksii/cbmdir/flow_drawing_tools/.config.sh'
alias flowcalcconfig='source /home/oleksii/cbmdir/flow_calculator/.config.sh'
alias qndiscrconfig='source /home/oleksii/cbmdir/qn_discriminator/.config.sh'

alias files644='find ./ -type f -print0 | xargs -0 chmod 644'
alias dirs755='find ./ -type d -print0 | xargs -0 chmod 755'

alias lx='pdflatex main'
alias bx='bibtex main'
alias latexclean='rm main.bbl main.aux main.log main.toc main.pdf main.out main.lot main.lof main.blg'
alias lbll='lx && bx && lx && lx'

alias rootl='root -l'
alias rootlq='root -l -q'
alias rootlb='root -l -b'
alias rootlbq='root -l -b -q'
alias r='root'
alias rl='rootl'
alias rlq='rootlq'
alias rlb='rootlb'
alias rlbq='rootlbq'

########################################################################################################

alias virgovae23='ssh -Y vae23.hpc.gsi.de'
alias virgo='ssh virgo.hpc.gsi.de'
alias virgo3='ssh -Y virgo3.hpc.gsi.de'
alias virgo0552='ssh lxbk0552'
alias virgo1130='ssh lxbk1130'
alias virgo1131='ssh lxbk1131'
alias virgo1132='ssh lxbk1132'
alias virgo1133='ssh lxbk1133'
alias singunohome='singularity shell -B /data.local1 -B /lustre --no-home /lustre/alice/singularity/singularity_o2dev_alma9.sif'
alias ssbs='cd; source .start-scratch-build-session.sh'
alias abbo2o='aliBuild build O2Physics'
alias aeo2p='cd /scratch/alice/lubynets/alice; alienv enter O2Physics/latest ninja/latest --shellrc'
alias sscratch='source /lustre/alice/containers/start-scratch-build-session.sh'
alias o2dep='$O2PHYSICS_ROOT/share/scripts/find_dependencies.py -t'
alias nano="/lustre/cbm/users/lubynets/soft/nano/install/bin/nano"
alias mntlustre='sshfs virgo-centos7.hpc.gsi.de:/lustre/cbm/users/lubynets/ /u/lubynets/Mount/Lustre'
alias cmake='/lustre/cbm/users/lubynets/soft/CMake/install_3.21/bin/cmake'
alias lucbm='cd /lustre/cbm/users/lubynets'
alias lu='cd /lustre/alice/users/lubynets'
alias scal='cd /scratch/alice/lubynets'
alias mji='make -j24 install'
alias mj='make -j24'
alias sq='squeue -u lubynets'
alias wsq='watch squeue -u lubynets'
alias sqsc='scancel -u lubynets'
alias cd1='cd ..'
alias cd2='cd ../..'
alias cd3='cd ../../..'
alias cd4='cd ../../../..'
alias cd5='cd ../../../../..'
alias cd6='cd ../../../../../..'
alias cd7='cd ../../../../../../..'
alias cd8='cd ../../../../../../../..'
alias cd9='cd ../../../../../../../../..'
alias gst='git status'
alias grv='git remote -v'
alias gpl='git pull'
alias glo='git log'
alias gdf='git diff'
alias gch='git checkout'
alias tm='tmux'
alias ta='tmux attach'
alias ta0='tmux attach -t 0'
alias ta1='tmux attach -t 1'
alias ta2='tmux attach -t 2'
alias rmlw='rm -rf log/ workdir/'
alias rmwl='rm -rf log/ workdir/'
alias qnanalysisconfig='source /lustre/cbm/users/lubynets/soft/QnAnalysis/install/bin/QnAnalysisConfig.sh'

alias checkspace='lfs quota -h -u $USER /lustre/cbm/users/$USER/'
alias countfiles='/lustre/cbm/users/lubynets/.counter_script.sh'
alias dirssize='du -sh -- *'

alias ROOT='source /lustre/alice/users/lubynets/soft/root/install_6.32_cpp17/bin/thisroot.sh'
alias rootl='root -l'
alias rootlq='root -l -q'
alias r='root'
alias rl='rootl'
alias rlq='rootlq'

sast () {
  if [ ! -d "/lustre" ]; then
    echo "You're outside Virgo & Lustre"
  else
    nodename=$(uname -n)
    if [[ $nodename == "lxbk1132" || $nodename == "lxbk1133" ]]; then
      echo "You're at Virgo vae23, not suited for alice software"
    else
      if [ ! -d "/.singularity.d" ]; then
      echo "You're outside singularity container"
      echo "You're outside ALICE environment"
      else
        echo "You're inside singularity container"
        AON=$(type -p alien-token-info)
        if [ -z $AON ]; then
          echo "You're outside ALICE environment"
        else
          echo "You're inside ALICE environment"
        fi
      fi
    fi
  fi
}

sqdet () {
  sacct -j $1 -o JobID,User,JobName,Partition,MaxVMSize,MaxRSS,State
}
