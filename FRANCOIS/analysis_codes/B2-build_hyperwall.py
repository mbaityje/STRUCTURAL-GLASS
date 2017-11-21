## build_hyperwall
## build the hyperwall using SVMs
## Reading the descriptors, output = wall+ offset 
## 
from __future__ import print_function

import sys
import numpy as np
import module_filenamesHandler 

import sklearn.svm
from sklearn.model_selection import cross_val_score

######  
## plots: 
import module_importPlotParams
import matplotlib.pyplot as plt 
plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.unicode']=True          
plt.rcParams.update(module_importPlotParams.fparams())
plt.ion()   ## optional
######


print('usage:')
print('run  build_hyperwall.py  [filename with "AVA=descriptors..."]   [mode=1,2]  ')

filename = sys.argv[1]
rootname = module_filenamesHandler.get_rootname(filename)
suffix   = module_filenamesHandler.get_suffix(filename)
mode     = int(sys.argv[2])

cv_cross = 2
cv_cross = int(sys.argv[3])


Natoms = int(float(module_filenamesHandler.filename_parser(filename[:-4], 'N')))
seed=4242
if "AVA=descriptors" in filename :
    X_all_descriptors_array = np.loadtxt(filename)
else:
    print("no way. Use correct file conttaing the ddescriptors (prodcued with B1).")
    raise SystemExit
Nsamples = len(X_all_descriptors_array[:,0])/2
print("found ", Nsamples, " samples (hard, soft)" )
#hard_descriptors = X_all_descriptors_array[0:Nsamples,:]
#soft_descriptors = X_all_descriptors_array[Nsamples: ,:]
y_correct_labels = np.zeros(2*Nsamples)-1
y_correct_labels[Nsamples:] += 2 

#Nrestricted=Nsamples
#X_all_descriptors_array[Nrestricted:2*Nrestricted] = X_all_descriptors_array[-Nrestricted:]
#X_all_descriptors_array = X_all_descriptors_array[:2*Nrestricted]
#y_correct_labels[Nrestricted:2*Nrestricted] = y_correct_labels[-Nrestricted:]
#y_correct_labels = y_correct_labels[:2*Nrestricted]


Ndescriptors = len(X_all_descriptors_array[0])
family = int(module_filenamesHandler.filename_parser(filename[:-4], 'fam'))
print("The family of descriptors is    family#",family,"   and each point (example) has   Ndescriptors=", Ndescriptors, "   descriptors")


##########################
### normalizing the descriptors, so they all come in equally importantly.
descriptor_averages = np.mean(X_all_descriptors_array,0)
descriptor_stdDevs  = np.std( X_all_descriptors_array,0)
X_all_descriptors_array -= descriptor_averages
X_all_descriptors_array /= descriptor_stdDevs
print("the averages and stdDev of descriptors are now: 0 and 1 ")
#print("or precicely (DEBUG):", np.mean(X_all_descriptors_array,0), np.std( X_all_descriptors_array,0) ) 
##########################


Cvalue = 1.0
verbose = True
probability_estimates = False ## (only usefull when performing  predict_log_proba
shrinking= True # False    ## setting it to True may spped up optimization.. but wouthout it I guess results are more reliable ? Also it is not sued in th eother LinearSVC method, I guess, so for comaprison it is better to do without it.
tolerance = 1e-3
####	problem_parameters_cache_size = 2000;
#	problem_parameters_eps = 1.0e-3;
#	problem_parameters_C = pow(2.0,C);
#	problem_parameters_gamma = pow(2.0,gamma);
#	cout << "C = " << problem_parameters_C << endl;

####	problem_parameters_nr_weight = 0;
####	problem_parameters_weight_label = NULL;
####	problem_parameters_weight = NULL;
####	problem_parameters_probability = 0;

####	problem_parameters_shrinking = 1;
####	problem_parameters_nu = 0.5;
####	problem_parameters_degree = 3;
####	problem_parameters_coef0 = 0;
####	problem_parameters_p = 0.1;



### perform cross-validation ##
if mode == 1:
    print("Mode 1 : Cross validation, i.e. finding the optimal \gamma, C, Nsamples values for the training set.")
    
    Kernel='linear'
#    clf = sklearn.svm.SVC(C=Cvalue, kernel=Kernel, probability=probability_estimates, shrinking=shrinking, tol=tolerance, cache_size=2000, verbose=verbose, max_iter=-1, decision_function_shape='ovr', random_state=seed)
    ### clf.support_vectors_    ## only supported using SVC, not LinearSVC ##

    clf = sklearn.svm.LinearSVC(C=Cvalue, penalty='l2', loss='squared_hinge', dual=False, tol=tolerance,  intercept_scaling=1, verbose=True, random_state=seed, multi_class='ovr', max_iter=1000000)

#    for Cvalue in [10**i for i in range(-4,3)]:
    for Cvalue in [100**i for i in range(-1,2)]: 
        clf.set_params(C=Cvalue)
        clf.set_params(verbose=False)


    #    Nrestricted = 1000
    #    X_restricted = X_all_descriptors_array[:2*Nrestricted]
    #    X_restricted[Nrestricted:] = X_all_descriptors_array[Nsamples:Nsamples+Nrestricted]
    #    y_resctricted= y_correct_labels[:2*Nrestricted]
    #    y_resctricted[Nrestricted:] = y_correct_labels[Nsamples:Nsamples+Nrestricted]

    #    clf.fit(X=X_restricted, y = y_resctricted)
    #    clf.fit(X=X_all_descriptors_array, y = y_correct_labels)
        ## dispaly support vectors:
        scores = sklearn.model_selection.cross_val_score(clf, X=X_all_descriptors_array, y=y_correct_labels , cv=cv_cross)
        
        print('C=', Cvalue, "  average=", str(np.mean(scores))[:4],  "  stdDev=", np.std(scores), " cv_cross=",cv_cross)#  "  scores=", scores,  







if mode ==2 : 
    Cvalue = 1
    clf = sklearn.svm.LinearSVC(C=Cvalue, penalty='l2', loss='squared_hinge', dual=False, tol=tolerance,  intercept_scaling=1, verbose=True, random_state=seed, multi_class='ovr', max_iter=1000000)
    clf.fit(X_all_descriptors_array, y_correct_labels)
    scores = clf.score(X_all_descriptors_array, y_correct_labels)
    print('C=', Cvalue, "  self score = ", scores , "   (cheating because using the training set as test set !!)")

    ## saving the model to a file: ##
    a_hyperplane = clf.coef_[0]
    b_intercept = clf.intercept_[0]
    hyperplane = np.zeros( (len(a_hyperplane), 3))
    hyperplane[:,0]= descriptor_averages
    hyperplane[:,1]= descriptor_stdDevs
    hyperplane[:,2]= a_hyperplane
    header= "saved SVM hyperplane and noramlization data.\nb = "+str(b_intercept)+" \nNOTE: the second line needs to stay in this format, in particular the last piece (space-separated) needs to be b itself. It is read using: b = float((flow.readline()).split()[-1]) \n mean(G)                  sigma(G)                 a_hyperplane (vector)"
    hyperplane_name = filename[:-4]+"_hyperplane.dat"
    np.savetxt(hyperplane_name, hyperplane, header=header)
    print("saved the hyperplane to file: \n", hyperplane_name)


#if mode == 3 :
#    
#    filename = filename[:-4]+"_hyperplane.dat"
#    
#    
#    hyperplane = np.loadtxt(filename)
#    descriptor_averages = hyperplane[:,0]
#    descriptor_stdDevs  = hyperplane[:,1]
#    a                   = hyperplane[:,2]
#    with open(filename, 'r') as flow:
#        flow.readline()
#        b = float((flow.readline()).split()[-1])

#    if "AVA=descriptors" in filename :
#        X_all_descriptors_array = np.loadtxt(filename)
#        X_all_descriptors_array -= descriptor_averages
#        X_all_descriptors_array /= descriptor_stdDevs
#    else:
#        print("no way. Use correct file conttaing the ddescriptors (prodcued with B1).")
#        raise SystemExit

#    ### DEBUG : check ######
#    y_calc = np.zeros(Nsamples*2)
#    for i in range(Nsamples*2):
#        y_calc[i] = np.dot(a,X_all_descriptors_array[i])+b

#    y_clf = clf.decision_function(X_all_descriptors_array)
#    print(y_clf - y_calc)
#    
#    hard_hist, bins = np.histogram(y_calc[:Nsamples], 100)
#    plt.plot(bins[:-1], hard_hist)
#    
#    soft_hist, bins = np.histogram(y_calc[Nsamples:], 100)
#    plt.plot(bins[:-1], soft_hist)
#    
#    total_hist, bins = np.histogram(y_calc[:], 100)
#    plt.plot(bins[:-1], total_hist)
#    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
    
#	
#	
#	
#	
#	
#import cPickle
## save the classifier
#with open('my_dumped_classifier.pkl', 'wb') as fid:
#    cPickle.dump(gnb, fid)    

## load it again
#with open('my_dumped_classifier.pkl', 'rb') as fid:
#    gnb_loaded = cPickle.load(fid)
#	


