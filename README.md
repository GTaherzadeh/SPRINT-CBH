# SPRINT-CBH
Carbohydrate-binding proteins play significant roles in many diseases including cancer. Here, we established a machine-learning-based method (called Sequencebased Prediction of Residue-level INTeraction sites of carbohydrates, SPRINT-CBH) to predict carbohydrate-binding sites in proteins by using Support Vector Machines (SVM). We found that integrating evolution-derived sequence profiles with additional information of sequence and predicted solvent accessible surface area leads to a reasonably accurate, robust, predictive method, with area under receiver operating characteristic curve (AUC) of 0.78 and 0.77, and Matthew’s correlation coefficient of 0.34 and 0.29, respectively for ten-fold cross validation and independent test without balancing binding and non-binding residues. The quality of the method is further demonstrated by having statistically significantly more binding residues predicted for carbohydrate-binding proteins than presumptive non-binding proteins in the human proteome, and by the bias of rare alleles toward predicted carbohydrate-binding sites for non-synonymous mutations from the 1000 genome project. 

Cite: Taherzadeh, G., Zhou, Y., Liew, A. W. C., & Yang, Y. (2016). Sequence-based prediction of protein–carbohydrate binding sites using support vector machines. Journal of chemical information and modeling, 56(10), 2115-2122.

Instruction:

* Protein-carbohydrate dataset are stored in Data.zip. Dataset file contains protein sequences labeled as 1 and 0 for binding and non-binding residues, respectively, and specified as train and test folders. Each pdb.txt contains actual binding residues used in this study. The pdb files are available in Data.zip.
* Run ./SPRINT-CBH.py for feature extraction and carbohydrate binding site prediction.
* The pre-trained model is available in Fullmodel.zip.
