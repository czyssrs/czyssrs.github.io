rm -r ours
cp -r ../runs/20170314_01/output ours
python compare_motif.py meme ours
