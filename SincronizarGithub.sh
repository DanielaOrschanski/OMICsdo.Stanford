git config --global user.name "DanielaOrschanski"
git config --global user.email "dorschanski@gmail.com"

sudo apt-get update

#token omicsdo
#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

cd /media/8tb01/Daniela/OMICsdo.Stanford
git init
git --version

#Cada vez que cambio algo
git add .
git commit -m "Actualizacion 30/9"
git push -u origin main

#si da error:
git pull origin main

#Si elimine cosas por espacio insuficiente:
sudo apt install git-filter-repo
git filter-repo --force \
  --path data/Pfam-A.hmm \
  --path data/Pfam-A.hmm.h3p \
  --path data/Pfam-A.hmm.h3f \
  --path data/Pfam-A.hmm.h3m \
  --invert-paths
git push origin main --force
git remote -v

#esto solo la primwra vez
git config user.name
git config user.email

git remote add origin https://github.com/DanielaOrschanski/OMICsdo.Stanford.git
git remote -v

git add .
git commit -m "Initial commit"
git branch -M main
git push -u origin main


##### ERROR
git clone --mirror https://github.com/DanielaOrschanski/OMICsdo.Stanford.git
cd OMICsdo.Stanford.git
sudo apt install git-filter-repo
#Guada14514213
git filter-repo --invert-paths --path SincronizarGithub.sh

#####
#tuve que eliminar archov grande porque me daba error:
rm Log.out
git filter-branch --force --index-filter \
'git rm --cached --ignore-unmatch Log.out' \
--prune-empty --tag-name-filter cat -- --all
#-----




################################################
#token omicsdo
#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


#esto solo la primwra vez
git add .
git commit -m "Initial"
git branch -M main
git pull --rebase origin main

git push -u origin main

