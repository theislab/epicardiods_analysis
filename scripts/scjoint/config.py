import torch
import os

class Config(object):
    def __init__(self):
        DB = '10x'
        self.use_cuda = True
        self.threads = 1

        if not self.use_cuda:
            self.device = torch.device('cpu')
        else:
            self.device = torch.device('cuda:0')
        
        if DB == '10x':
            # DB info
            self.number_of_class = 31
            self.input_size = 16949
            data_path = '/storage/groups/ml01/workspace/laura.martens/moretti_colab/scjoint'
            self.rna_paths = [os.path.join(data_path, 'exprs_10xPBMC_rna.npz')]
            self.rna_labels = [os.path.join(data_path, 'cellType_10xPBMC_rna.txt')]		
            self.atac_paths = [os.path.join(data_path, 'exprs_10xPBMC_atac.npz')]
            self.atac_labels = [] #Optional. If atac_labels are provided, accuracy after knn would be provided.
            self.rna_protein_paths = []
            self.atac_protein_paths = []
            
            # Training config            
            self.batch_size = 256
            self.lr_stage1 = 0.01
            self.lr_stage3 = 0.01
            self.lr_decay_epoch = 20
            self.epochs_stage1 = 20
            self.epochs_stage3 = 20
            self.p = 0.8
            self.embedding_size = 64
            self.momentum = 0.9
            self.center_weight = 1
            self.with_crossentorpy = True
            self.seed = 1
            self.checkpoint = ''


            

        



