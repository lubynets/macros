alias gsi='ssh lxpool'
# alias mntgsi='sshfs lubynets@lxi001.gsi.de:/u/lubynets/ /home/oleksii/Mount/gsi/'
alias mntgsi="sshfs -o ssh_command='ssh -J lubynets@lxlogin.gsi.de' lubynets@lxpool.gsi.de:/u/lubynets/ /home/oleksii/Mount/gsi/"
alias mntlustre2steps='sshfs lubynets@lxi001.gsi.de:/u/lubynets/Mount/Lustre/ /home/oleksii/Mount/Lustre/'
alias mntlustre="sshfs -o ssh_command='ssh -J lubynets@lxi001.gsi.de' lubynets@lustre.hpc.gsi.de:/lustre/cbm/users/lubynets/ /home/oleksii/Mount/Lustre/"
alias mntlustrealice="sshfs -o ssh_command='ssh -J lxpool' lubynets@virgo3.hpc.gsi.de:/lustre/alice/users/lubynets/ /home/oleksii/Mount/Lustre_alice/"
alias mntlustrecbm="sshfs -o ssh_command='ssh -J lxpool' lubynets@virgo3.hpc.gsi.de:/lustre/cbm/users/lubynets/ /home/oleksii/Mount/Lustre/"
alias mntscratchalice="sshfs -o ssh_command='ssh -J lxpool' lubynets@virgo.hpc.gsi.de:/scratch/alice/lubynets/ /home/oleksii/Mount/scratch_alice/"
alias mmg="mntgsi && mntlustrealice && mntscratchalice && gsi"
alias mntall='mntgsi && mntlustrealice && mntscratchalice'
alias umntgsi='sudo umount ~/Mount/gsi'
alias umntlustre='sudo umount ~/Mount/Lustre_alice'
alias umntlustrealice='sudo umount ~/Mount/Lustre'
alias umntscratchalice='sudo umount ~/Mount/scratch_alice'
alias umntall='umntgsi; umntlustrealice; umntscratchalice'
alias kssh='killall -9 ssh'
alias ws='watch sensors'

alias cern='ssh -X ollubyne@lxplus981.cern.ch'
alias mntcern='sshfs ollubyne@lxplus.cern.ch:/afs/cern.ch/user/o/ollubyne/ /home/oleksii/Mount/cern/'

alias mji='make -j3 install'
alias mj='make -j3'
alias mi='make install'
alias m='make'
alias mjit='mji; make test'
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
alias gadd='git add .'
alias gcm='git commit -m'
alias kt='killall -9 telegram'
alias my-clang-format='~/.clang-format'
alias pdf2png='pdftoppm -png -cropbox -singlefile'
alias pdf2pngall='pdftoppm -png -cropbox'
pdf2pngpage() {
  pdftoppm -png -cropbox -singlefile -l $1 -f $1 $2 $3
}
alias cpdf='~/.cpdf' # removes pictures from pdf: cpdf -draft input.pdf -o output.pdf
MPAD() {
  if [ -n "$2" ]; then
  EXE="-o $2"
  fi
  echo "g++ $1 ~/cbmdir/macros_on_git/MultiPad/MultiPicture.cpp -I ~/cbmdir/macros_on_git/MultiPad $EXE"
  g++ $1 ~/cbmdir/macros_on_git/MultiPad/MultiPicture.cpp -I ~/cbmdir/macros_on_git/MultiPad $EXE
}
rena() {
  rename "s/$1/$2/" *
}

alias ROOT24='source /home/oleksii/soft/Root/install_root_6.24/bin/thisroot.sh'
alias ROOT32='source /home/oleksii/soft/Root/install_root_6.32/bin/thisroot.sh'
alias AT='source /home/oleksii/soft/AnalysisTree/install_230/bin/AnalysisTreeConfig.sh'
alias qntoolsconfig='source /home/oleksii/soft/.qntoolsconfig.sh'
alias flowdrawconfig='source /home/oleksii/cbmdir/flow_drawing_tools/.config.sh'
alias flowcalcconfig='source /home/oleksii/cbmdir/flow_calculator/.config.sh'
alias qndiscrconfig='source /home/oleksii/cbmdir/qn_discriminator/.config.sh'
alias qa2config='source /home/oleksii/alidir/macros_on_git/qa2/exe_based/.config.sh'

alias cmakedebug='cmake -DCMAKE_BUILD_TYPE=DEBUG ..'
alias cmakerelease='cmake -DCMAKE_BUILD_TYPE=RELEASE ..'

alias files644='find ./ -type f -print0 | xargs -0 chmod 644'
alias dirs755='find ./ -type d -print0 | xargs -0 chmod 755'
alias dirssize='du -sh -- *'

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
alias virgovae25='ssh -Y vae25.hpc.gsi.de'
alias virgovae='virgovae25'
alias virgo='ssh virgo.hpc.gsi.de'
alias virgo3='ssh -Y virgo3.hpc.gsi.de'
alias virgo1130='ssh lxbk1130'
alias virgo1131='ssh lxbk1131'
alias virgo1132='ssh lxbk1132'
alias virgo1133='ssh lxbk1133'
alias virgo30='virgo1130'
alias virgo31='virgo1131'
alias virgo32='virgo1132'
alias virgo33='virgo1133'
alias virgo36='ssh lxbk1136'
alias virgo37='ssh lxbk1137'
alias bdtsif='srun --pty --reservation="cluster_test" --partition=main --mem=46G --cpus-per-task=24 -t 08:00:00 -- apptainer shell -B /u /lustre/alice/users/lubynets/singularities/bdt.sif'
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
alias mi='make install'
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
alias gcm='git commit -m'
alias tm='tmux'
alias ta='tmux attach'
alias ta0='tmux attach -t 0'
alias ta1='tmux attach -t 1'
alias ta2='tmux attach -t 2'
alias ktm='tmux kill-session'
alias rmlw='rm -rf log/ workdir/'
alias rmwl='rm -rf log/ workdir/'
alias ee='exit'

alias checkspace='lfs quota -h -u $USER /lustre/alice/users/$USER/'
alias countfiles='/lustre/cbm/users/lubynets/.counter_script.sh'
alias dirssize='du -sh -- *'

alias cmakedebug='cmake -DCMAKE_BUILD_TYPE=DEBUG ..'
alias cmakerelease='cmake -DCMAKE_BUILD_TYPE=RELEASE ..'

alias ROOT='source /lustre/alice/users/lubynets/soft/root/install_6.32_cpp17_vae25/bin/thisroot.sh'
alias ROOT_vae23='source /lustre/alice/users/lubynets/soft/root/install_6.32_cpp17/bin/thisroot.sh'
alias rootl='root -l'
alias rootlq='root -l -q'
alias r='root'
alias rl='rootl'
alias rlq='rootlq'

alias tarlist='tar -ztvf'
alias tarextract='tar -zxvf'

sast () {
  if [ ! -d "/lustre" ]; then
    echo "You're outside Virgo & Lustre"
  else
    nodename=$(uname -n)
    if [[ $nodename == "lxbk1132" || $nodename == "lxbk1133" || $nodename == "lxbk1136" ||$nodename == "lxbk1137" ]]; then
      echo "You're at Virgo vae23 or vae25, not suited for alice software"
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
