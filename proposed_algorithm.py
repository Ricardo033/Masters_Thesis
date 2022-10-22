#-----------------------------------------------------#
# Propposed method for unsupervised sample selection  #
# by: Ricardo Andres Castañeda Rueda                  #
# student id: r0731529                                #
#-----------------------------------------------------#

#====================#
# Required packages  #
#====================#

import pandas as pd
import numpy as np
import random as rnd
import itertools
import numpy as geek
import matplotlib.pyplot as plt


class domain_invariant_unsupervised_sample_selection(object):

# ---------------------------- Z matrix with subset of samples from X --------------------------------        
    def sub_matrix(x_matrix, num_samples):
            ''' 
This function selects rows or samples from a population matrix (X) randomly. The goal is to obtain a submatrix matrix (Z) that contains a subset of samples to use in a calibration model. Besides, this function yields another matrix (Y) that contains the samples that were not included in the submatrix matrix (Z) and a tow vectors (Z_samples, Y_samples) that indicate the samples’ index contained in the resulted matrices.    
        
        --- Input ---
        
        x_matrix: Population matrix that includes all samples to evaluate as possible candidates to use in the calibration model.      
        num_samples: The number of samples the user wants to include in the calibration model. This variable indicates how many samples are to be included in the submatrix (Z).   
                
        --- Output ---
        Z: An array that contains some samples from the population matrix (X). The length of this matrix (Z) depends on the variable num_samples, and it's expected to be smaller than (X), but it can be set depending on the user's needs. 
        Y: An array that contains the samples from the population matrix (X), that where not included in the submatrix (Z).
        Z_samples: A one-dimensional array that shows the index of the samples contained in the submatrix (Z). Please note that this index shows the position that a sample has in the population matrix (X), and this way, it is possible to identify which samples were included in the submatrix (Z). 
         Y_samples: A one-dimensional array that shows the index of the samples not contained in the submatrix (Z). Then it identifies which samples were not included in the submatrix (Z). 
        
        '''
        
            Z_samples = np.zeros(num_samples).astype(np.int32)
            Z= np.zeros((num_samples, len(x_matrix[1])))
            Y= x_matrix.copy()
            Num_rows = num_samples
            i = 0
            l= len(x_matrix)
            index = np.array(range(0,l))
            while i < Num_rows:
                r = rnd.randint(0,l-1)
                if r not in Z_samples: 
                    Z_samples[i] = r  
                    row = np.copy(x_matrix[r,:])
                    Z[i] =  row
                    i += 1

            Y_samples = np.delete(index.copy(), Z_samples, axis = 0)  
            Y = np.delete(x_matrix.copy(), Z_samples, axis = 0)

            result=[Z,Y,Z_samples, Y_samples]

            print("Z= element [0]", "dimetion:", result[0].shape,"numpy.ndarray")
            print("Y= element [1]", "dimetion:", result[1].shape,"numpy.ndarray")
            print("Zsamples= element [2]", "length:",len(result[2]),"numpy.ndarray")
            print("Ysamples= element [3]", "length:",len(result[3]),"numpy.ndarray")
            return result

 #------
    def diuss_max(X,Z,Y,Z_samples,Y_samples,iterations):

        #Mathematical criterion
        def criterion(x_matrix, z_matrix):

            cov_X = np.cov(x_matrix.T,ddof=0)
            cov_Z = np.cov(z_matrix.T,ddof=0)

            D = (cov_X-cov_Z)
            svd_D = np.linalg.svd(D)
            crt = max(abs(svd_D[1]))
            return crt
        # Vector of samples in Binary form
        def binary (x_matrix, z_samples):

            index = np.array(range(0,len(x_matrix)))
            samples = z_samples
            vector = np.zeros(len(index))
            for i in range(len(samples)):
                for j in range(len(index)):

                    if samples[i] == index[j]:
                        vector[j] = 1
            return(vector) 

        #Matrices
        Z_new = np.zeros(( len(Z), len(Z[1])))
        Y_new = np.zeros(( len(Y), len(Y[1])))
        Z_prev = Z.copy()
        Y_prev = Y.copy()

        #List of samples contained in the matrices
        Z_Samples_new = np.zeros((0,len(Z))).astype(np.int32) #**
        Y_Samples_new = np.zeros((0,len(Y))).astype(np.int32) #**
        Z_Samples_prev = Z_samples.copy() #**
        Y_Samples_prev = Y_samples.copy() #**

        # Number of samples or rows in the matrices
        rows_X = len(X)
        rows_Y = len(Y)
        rows_Z = len(Z)
        # Cov matrices and mathematical criterion
        crit = criterion(x_matrix=X, z_matrix=Z)
        all_crit = np.zeros(iterations+1)
        all_crit[0] = crit

        for i in range(iterations):
            # Random rows
            r1 = rnd.randint(0,rows_Z-1) 
            r2 = rnd.randint(0,rows_Y-1) 
            # Swaping variables in Z
            Z_new = np.delete(Z_prev.copy(), r1, axis = 0)
            Z_new = np.insert(Z_new, r1, Y_prev[r2], axis =0)

            ## Swaping variables in Y
            Y_new = np.delete(Y_prev.copy(), r2, axis = 0)
            Y_new = np.insert(Y_new, r2, Z_prev[r1], axis =0)

            #Step 2: Criterion
            New_crit = criterion(x_matrix=X, z_matrix=Z_new)

            #Step 3: Selecting Best matrix
            if New_crit < crit:
                crit = New_crit
                Z_prev = Z_new.copy()
                Y_prev = Y_new.copy()

                # List of samples in the matrices
                Z_Samples_new = np.delete(Z_Samples_prev.copy(), r1, axis = 0)
                Z_Samples_new = np.insert(Z_Samples_new, r1, Y_Samples_prev[r2], axis =0)
                Y_Samples_new = np.delete(Y_Samples_prev.copy(), r2, axis = 0)
                Y_Samples_new = np.insert(Y_Samples_new, r2, Z_Samples_prev[r1], axis =0)
                Z_Samples_prev = Z_Samples_new.copy() 
                Y_Samples_prev = Y_Samples_new.copy()

                # Vector containing all calculated criterions
            all_crit[i+1] = New_crit

        samples_binary = binary(x_matrix=X, z_samples=Z_Samples_prev)  
        samples_binary = samples_binary.astype(np.int32)
        result = [Z_prev, Z_Samples_prev, Y_prev, Y_Samples_prev, all_crit, samples_binary]
        return result
 #--------------
    def diuss_sum(X,Z,Y,Z_samples,Y_samples,iterations):

        #Mathematical criterion
        def criterion(x_matrix, z_matrix):

            cov_X = np.cov(x_matrix.T,ddof=0)
            cov_Z = np.cov(z_matrix.T,ddof=0)

            D = (cov_X-cov_Z)
            svd_D = np.linalg.svd(D)
            crt2 = sum(abs(svd_D[1]))
            return crt2
        # Vector of samples in Binary form
        def binary (x_matrix, z_samples):

            index = np.array(range(0,len(x_matrix)))
            samples = z_samples
            vector = np.zeros(len(index))
            for i in range(len(samples)):
                for j in range(len(index)):

                    if samples[i] == index[j]:
                        vector[j] = 1
            return(vector) 

        #Matrices
        Z_new = np.zeros(( len(Z), len(Z[1])))
        Y_new = np.zeros(( len(Y), len(Y[1])))
        Z_prev = Z.copy()
        Y_prev = Y.copy()

        #List of samples contained in the matrices
        Z_Samples_new = np.zeros((0,len(Z))).astype(np.int32) #**
        Y_Samples_new = np.zeros((0,len(Y))).astype(np.int32) #**
        Z_Samples_prev = Z_samples.copy() #**
        Y_Samples_prev = Y_samples.copy() #**

        # Number of samples or rows in the matrices
        rows_X = len(X)
        rows_Y = len(Y)
        rows_Z = len(Z)
        # Cov matrices and mathematical criterion
        crit = criterion(x_matrix=X, z_matrix=Z)
        all_crit = np.zeros(iterations+1)
        all_crit[0] = crit

        for i in range(iterations):
            # Random rows
            r1 = rnd.randint(0,rows_Z-1) 
            r2 = rnd.randint(0,rows_Y-1) 
            # Swaping variables in Z
            Z_new = np.delete(Z_prev.copy(), r1, axis = 0)
            Z_new = np.insert(Z_new, r1, Y_prev[r2], axis =0)

            ## Swaping variables in Y
            Y_new = np.delete(Y_prev.copy(), r2, axis = 0)
            Y_new = np.insert(Y_new, r2, Z_prev[r1], axis =0)

            #Step 2: Criterion
            New_crit = criterion(x_matrix=X, z_matrix=Z_new)

            #Step 3: Selecting Best matrix
            if New_crit < crit:
                crit = New_crit
                Z_prev = Z_new.copy()
                Y_prev = Y_new.copy()

                # List of samples in the matrices
                Z_Samples_new = np.delete(Z_Samples_prev.copy(), r1, axis = 0)
                Z_Samples_new = np.insert(Z_Samples_new, r1, Y_Samples_prev[r2], axis =0)
                Y_Samples_new = np.delete(Y_Samples_prev.copy(), r2, axis = 0)
                Y_Samples_new = np.insert(Y_Samples_new, r2, Z_Samples_prev[r1], axis =0)
                Z_Samples_prev = Z_Samples_new.copy() 
                Y_Samples_prev = Y_Samples_new.copy()

                # Vector containing all calculated criterions
            all_crit[i+1] = New_crit

        samples_binary = binary(x_matrix=X, z_samples=Z_Samples_prev)  
        samples_binary = samples_binary.astype(np.int32)
        result = [Z_prev, Z_Samples_prev, Y_prev, Y_Samples_prev, all_crit, samples_binary]
        return result
 #---------------------------------- Count of ones for vectors of selected samples --------------
    def count_of_ones(array):
        a=0
        for i in range(len(array)): 
            if  array[i] == 1:
                a = a +1
        return(a)
#-----------------
# Function to plot a vector that contain the recorde result from the criterion obtained at each iteration
    def crit_behavior(all_crit):
        a = all_crit
        index=np.array(range(0,len(a)))
        plt.plot(index,a)
        plt.ylabel('All sorted criteria')
        plt.show()
        return(plt.show())
