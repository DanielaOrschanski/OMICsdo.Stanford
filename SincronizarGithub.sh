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
git remote add origin https://github.com/DanielaOrschanski/OMICsdo.Stanford.git
git remote -v
#git remote remove origin
git status
git branch
git push -u origin main

#tuve que eliminar archov grande porque me daba error:
rm Log.out
git filter-branch --force --index-filter \
'git rm --cached --ignore-unmatch Log.out' \
--prune-empty --tag-name-filter cat -- --all
#-----




################################################
#token omicsdo
#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

cd /media/4tb2/Daniela/Biota/PipelineBiota
git init
git --version

git remote add origin https://github.com/DanielaOrschanski/PipelineBiota.git
git remote -v
git remote remove origin


#esto solo la primwra vez
git add .
git commit -m "Initial"
git branch -M main
git pull --rebase origin main

git push -u origin main

