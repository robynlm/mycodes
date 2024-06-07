import numpy as np
import math
import scipy.io as sio
import matplotlib.pyplot as plt
import copy
plt.close('all')

#to do next make the initialise and forward functions

class CNN:
    """"""
    def __init__(self, filter_structure, pool_structure, neural_structure, image_shape):
        """"""
        self.nbr_jumps = 1
        self.image_shape = [image_shape]
        self.filter_structure = filter_structure
        self.filter_weight = []
        self.filter_bias = []
        # nbr filters, size, stride, padding
        self.pool_structure = pool_structure
        # nbr pooling, size, stride, function
        self.neural_structure = neural_structure
        self.neural_weight = []
        self.neural_bias =[]
        
    ##########################################################################
    #                              Initialise
    ##########################################################################
        
    def initialise_network(self, train_data):
        """ Initialise the network """
        conv_output = self.initial_conv(train_data)
        features_vector = [np.ravel(conv_output[i])
                           for i in range(np.shape(conv_output)[0])]
        print('   Filter part done')
        
        self.neural_structure = np.append(np.shape(features_vector[0]),
                                          self.neural_structure)
        for i in range(len(self.neural_structure)-1):
            self.neural_weight += [np.random.randn(self.neural_structure[i+1],
                                   self.neural_structure[i])]
            self.neural_bias += [np.random.randn(self.neural_structure[i+1])]
        print('   Neural part done')

    def initial_conv(self, images):
        """ Initialise the filters """
        filter_weighted_input = []
        nbr_layers = len(self.filter_structure)
        for i in range(nbr_layers):
            filter_layer = self.filter_structure[i]
            feature_filter = []
            for j in range(filter_layer['nbr filters']):
                feature_filter += [np.random.randn(np.shape(images[0])[0],
                                   filter_layer['size'], filter_layer['size'])]
            self.filter_weight += [feature_filter]
            self.filter_bias += [[np.random.rand()
                                  for k in range(filter_layer['nbr filters'])]]
            
            images = self.filterconv(images, filter_layer,
                                     self.filter_weight[i], self.filter_bias[i])
            
            filter_weighted_input += [images]
            images = [self.factive(np.array(images[a])) for a in range(len(images))]
            self.image_shape += [np.shape(images[0])]
            pool_layer = self.pool_structure[i]
            for j in range(pool_layer['nbr pooling']):
                images = [self.pool(images[m], pool_layer) for m in range(len(images))]
                self.image_shape += [np.shape(images[0])]
        return images
        
    ##########################################################################
    #                                 Forward
    ##########################################################################
        
    def forward(self, data):
        filter_outputs, filter_weighted_input = self.conv(data)
        features = filter_outputs[2]
        features_vector = [np.ravel(features[i]) for i in range(np.shape(features)[0])]
        aL, neural_output, neural_weighted_input = self.network(features_vector)
        return (filter_outputs, filter_weighted_input, features_vector,
                aL, neural_output, neural_weighted_input)

    ##########################################################################
    #                        Image extraction part
    ##########################################################################

    def padding(self, images, filter_padding):
        n = filter_padding
        padd2 = filter_padding*2
        Len = np.shape(images)
        C = n+Len[2]
        D = n+Len[3]
        padded = np.zeros([Len[0], Len[1], Len[2]+padd2, Len[3]+padd2])
        for i in range(len(images)):
            padded[i, :, n:C, n:D] = images[i]
        return padded

    def convolution(self, padded_images, layer, filters, bias):
        Len = np.shape(padded_images)
        size = layer['size']
        stride = layer['stride']
        n = (Len[2]-size)/float(stride)+1

        if n != int(n):
            print('ERROR convolution output size = ', n)
            print('')
        n = int(n)
        outcome = []
        for k in range(np.shape(filters)[0]):
            Filter = filters[k]
            outcomek = np.zeros([n, n])
            for i in range(n):
                di = stride*i
                maxi = di+size
                for j in range(n):
                    dj = stride*j
                    maxj = dj+size
                    outcomek[i,j] = sum(sum(sum(padded_images[:, di:maxi, dj:maxj]
                                                * Filter))) + bias[k]
            outcome += [outcomek]
        return outcome

    def filterconv(self, images, filter_layer, filters, bias):
        padded_images = self.padding(images, filter_layer['padding'])
        convolved = []
        for i in range(len(padded_images)):
            convolved += [self.convolution(padded_images[i], filter_layer, filters, bias)]
        return convolved

    def pool(self, income, pool_layer):
        function = pool_layer['function']
        if function == 'max':
            def func(A):
                return np.amax(A)
        elif function == 'mean':
            def func(A):
                return np.mean(A)
        else:
            print('Function not well defined')
        size = int(pool_layer['size'])
        stride = int(pool_layer['stride'])
        outcome_size = int(round((len(income[0])-size)/float(stride)+1))
        outcome = np.zeros([np.shape(income)[0], outcome_size, outcome_size])
        for k in range(np.shape(income)[0]):
            for i in range(outcome_size):
                di = stride*i
                for j in range(outcome_size):
                    dj = stride*j
                    matrix = income[k][di:di+size, dj:dj+size]
                    outcome[k][i][j] = func(matrix)
        return outcome

    def conv(self, images):
        output = [images]
        weighted_input = []
        nbr_layers = len(self.filter_structure)
        for i in range(nbr_layers):
            filter_layer = self.filter_structure[i]
            images = self.filterconv(images, filter_layer,
                                     self.filter_weight[i], self.filter_bias[i])
            weighted_input += [images]
            images = [self.factive(np.array(images[k])) for k in range(len(images))]
            pool_layer = self.pool_structure[i]
            for j in range(int(pool_layer['nbr pooling'])):
                images = [self.pool(images[k], pool_layer) for k in range(len(images))]
            output += [images]
        return output, weighted_input

    ##########################################################################
    ########################  Neural network part  ###########################
    ##########################################################################
        
    def neural_layer_output(self, w, x, b):
        z = np.array(w).dot(x) + np.array([b]).T
        a = self.active(np.array(z))
        return a, z

    def active(self, z):
        z = np.clip( z, -500, 500 )
        a = 1./(1.+np.exp(-z))
        return a

    def dactive(self, z):
        a = self.active(z)
        r = a*(1-a)
        return r

    def factive(self, z):
        return np.maximum(0, z)

    def dfactive(self, z):
        s = np.sign(z)
        ex = (s+1)/2
        return ex

    def softmax(self, z):
        z_exp = [np.exp(i) for i in z]
        softmax = [i / sum(z_exp) for i in z_exp]
        return softmax

    def network(self, data):
        neural_output = []
        a_output = [[[] for i in range(len(self.neural_structure)-1)]
                    for i in range(len(data))]
        weighted_input = [[[] for i in range(len(self.neural_structure)-1)]
                          for i in range(len(data))]
        for i in range(len(data)):
            x = np.array([data[i]]).T
            for j in range(len(self.neural_structure)-1):
                dummy = self.neural_layer_output(self.neural_weight[j], x,
                                                 self.neural_bias[j])
                x = dummy[0]
                a_output[i][j] = np.ravel(x)
                weighted_input[i][j] = np.ravel(dummy[1])
            neural_output += [np.ravel(x)]
        return neural_output, a_output, weighted_input


    ##########################################################################
    #####################  Learning - Backpropagation  #######################
    ##########################################################################

    def cost(self, output, label, size):
        C = sum(sum(abs(np.array(output)-np.array(label))**2))/(2.*size)
        return C

    def upsample(self, A, Bsize, pool_layer):
        B = np.zeros([Bsize, Bsize])
        psize = float(pool_layer['size'])
        div = psize**2
        for i in range(Bsize):
            I = int(round(i/psize))
            for j in range(Bsize):
                B[i, j] += A[I, int(round(j/psize))]/div
        return B

    def unravel(self, B, dim):
        C = np.zeros(dim)
        for k in range(dim[0]):
            K = k*dim[1]*dim[2]
            for i in range(dim[1]):
                I = i*dim[2]
                for j in range(dim[2]):
                    C[k][i][j] = B[K+I+j]
        return C

    def kernel_error(self, activ, err, padd):
        err = self.padding(err, padd)
        errk = [[[] for q in range(np.shape(activ)[1])] for p in range(np.shape(err)[1])]
        for q in range(np.shape(err)[1]): #nbr filters on that layer
            for p in range(np.shape(activ)[1]): #depth
                dummy = np.zeros([fsize, fsize])
                for u in range(fsize):
                    for v in range(fsize):
                        vect1 = err[:, q, :, :]
                        vect1[:, :u, :v] = np.zeros(np.shape(vect1[:, :u, :v]))
                        vect2 = activ[:, p, :, :]
                        dummy[u,v] = np.einsum('kij, kij', vect1, vect2)
                errk[q][p] = dummy
        return errk

    def kernel_backprop_error(self, err, kn, paddval):
        errS = [[[] for p in range(np.shape(kn)[1])] for k in range(np.shape(err)[0])]
        lenij = np.shape(kn)[0]
        err = self.padding(err, paddval)
        for k in range(np.shape(err)[0]):
            for p in range(np.shape(kn)[1]):
                dummy = np.zeros([lenij, lenij])
                for i in range(lenij):
                    for j in range(lenij):
                        uv = np.shape(kn)[2]
                        vect1 = err[k, :, i:i+uv, j:j+uv]
                        vect2 = kn[:, p, :, :]
                        dummy[i][j] = np.einsum('quv, quv', vect1, vect2)
                errS[k][p] = dummy
        return errS


    def accuracy(self, classification, testerY):
        count = 0
        for i in range(len(testerY)):
            if np.argmax(classification[i]) == np.argmax(testerY[i]):
                count +=1
        return count/float(len(testerY))
        
    def jump(self):
        next_jump = input("Jump? (y/n) ")
        if next_jump == 'y':
            self.nbr_jumps += 1
            locnoise = 0.0
            scalenoise = float(input("What variance for the noise? (float) "))
            nbr_previous_epoch = 0
            
            for i in range(2):
                self.filter_weight[i] += np.random.normal(loc=locnoise,
                                                         scale=scalenoise,
                                                         size=np.shape(self.filter_weight[i]))
                self.filter_bias[i] += np.random.normal(loc=locnoise,
                                                       scale=scalenoise,
                                                       size=np.shape(self.filter_bias[i]))
                self.neural_weight[i] += np.random.normal(loc=locnoise,
                                                         scale=scalenoise,
                                                         size=np.shape(self.neural_weight[i]))
                self.neural_bias[i] += np.random.normal(loc=locnoise,
                                                       scale=scalenoise,
                                                       size=np.shape(self.neural_bias[i]))

##########################################################################
#                           Data
##########################################################################

N = 28
zeros = [[np.zeros(N) for i in range(N)] for j in range(4)]
zeros_label = np.array([1, 0])

ones = [[np.ones(N) for i in range(N)] for j in range(4)]
ones_label = np.array([0, 1])

nbr_images_train = 2
trainerX = [zeros, ones] * nbr_images_train
trainerY = [zeros_label, ones_label] * nbr_images_train
nbr_images_test = 1
testerX = [zeros, ones] * nbr_images_test
testerY = [zeros_label, ones_label] * nbr_images_test

##########################################################################
#                           Structure
##########################################################################

np.random.seed(2)
# ---------------- Filter
fsize = 5
layer1 = {'nbr filters': 6, 'size': fsize, 'stride': 1, 'padding': 0}
layer2 = {'nbr filters': 12, 'size': fsize, 'stride': 1, 'padding': 0}
filter_structure = [layer1, layer2]
 # nbr filters, size, stride, padding
backprop_padding = np.array([2, 2])

psize = 3
layer1 = {'nbr pooling': 1, 'size': 3, 'stride': 2, 'function': 'max'}
layer2 = {'nbr pooling': 1, 'size': 3, 'stride': 2, 'function': 'max'}
pool_structure = [layer1, layer2]
 # nbr pooling, size, stride, function

# ---------------- Network
neural_structure = np.array([64, len(trainerY[0])])
#not including input neurons (features) and keeping size of label as nbr of ouput neurons
learning_rate = 3e-1
learning_rate_decay = 0. #5e-4
nbr_epoch = 2

# ---------------- Baches
bach_size = 128
bach_size_ini = 8
nbr_bach_train = int(math.ceil(len(trainerY)/float(bach_size)))
nbr_bach_test = int(math.ceil(len(testerY)/float(bach_size)))

# ---------------- Weights and biases
WandB = 'initialise'
#WandB = 'previous'

# don't be stupid, stupid
CnCdiff = 0.

image_shape = np.shape(trainerX[0])
cnn = CNN(filter_structure, pool_structure, neural_structure, image_shape)
###############################################################################
#                        INITIAL WEIGHTS AND BIASES
###############################################################################
print('')
print('--------------------INITIAL WEIGHTS AND BIASES-------------------')

if WandB == 'initialise':
    epoch_nbr = 0
    nbr_previous_epoch = 0
    train_data = trainerX[0:bach_size_ini]
    train_label = trainerY[0:bach_size_ini]
    cnn.initialise_network(train_data)

if WandB == 'previous':
    info = sio.loadmat('WandB3e-4.mat')
    epoch_nbr = info['epoch'][0][0]
    cnn.filter_weights = info['filter_weights'][0]
    FBiases = info['filter_biases'][0]
    filter_biases = []
    for i in range(len(FBiases)):
        filter_biases += [FBiases[i][0]]
    cnn.image_shape = info['image_shape']
    cnn.neural_structure = info['neural_structure'][0]
    cnn.neural_weights = info['neural_weights'][0]
    NBiases = info['neural_biases'][0]
    neural_biases = []
    nbr_previous_epoch = 2
    for i in range(len(NBiases)):
        neural_biases += [NBiases[i][0]]

###############################################################################
#                                TRAINING
###############################################################################
baches_train = np.arange(nbr_bach_train)
bb = baches_train*bach_size
print('')
print('----------------------------TRAINING-----------------------------')

Costs = []
Crecord = []
step_nbr_record = []
J = 0
Acc = []
while J<cnn.nbr_jumps:
    C = 15.
    Cn = 10.
    epoch = 0.
    Bforbreak = 0
    t = nbr_bach_train*nbr_previous_epoch
    dl = learning_rate
    Nepoch = nbr_epoch
    while epoch<Nepoch:
        b = 0
        while b<nbr_bach_train and abs(Cn-C)>CnCdiff:
            train_data = trainerX[bb[b]:bb[b]+bach_size]
            train_label = trainerY[bb[b]:bb[b]+bach_size]
            filter_outputs, filter_weighted_input, features_vector, aL, neural_output, neural_weighted_input = cnn.forward(train_data)
            
            nbr_trainer = len(aL)
            al = [neural_output[i][0] for i in range(nbr_trainer)]
            zL = [np.ravel(neural_weighted_input[i][1]) for i in range(nbr_trainer)]
            zl = [np.ravel(neural_weighted_input[i][0]) for i in range(nbr_trainer)]
            
            C = Cn
            Cn = cnn.cost(aL, train_label, nbr_trainer)
            Crecord += [Cn]
            acc = cnn.accuracy(aL, train_label)
            Acc += [acc]
            print('\r' '    >  epoch =', int(epoch+1), '/', Nepoch,
                  ',  bach =', b+1, '/', nbr_bach_train,
                  ',  loss =', "%.5f" % np.mean(Crecord[sum(step_nbr_record):]),
                  ',  acc =', "%.4f" % np.mean(Acc[sum(step_nbr_record):]))

            # ------------------------------------------
            # ---------- Error in neural part ----------
            # ------------------------------------------
            # ouput layer
            errL = (aL - np.array(train_label)) * cnn.dactive(np.array(zL))
            if np.amax(np.abs(errL)) == 0.0:
                print('L, Gradient Lost')
                Bforbreak = 1
                break
            # hidden layer
            errl = (errL.dot(np.array(cnn.neural_weight[1]))) * cnn.dactive(np.array(zl))
            if np.amax(np.abs(errl)) == 0.0:
                print('l, Gradient Lost')
                Bforbreak = 1
                break
            # 2nd convolutional layer
            errf2 = (errl.dot(np.array(cnn.neural_weight[0])))
                    # go back to tensor format
            errF2 = [cnn.unravel(errf2[i], cnn.image_shape[4])
                     for i in range(nbr_trainer)]
                    # unpool
            errc2 = [[cnn.upsample(errF2[j][i], cnn.image_shape[3][1], pool_structure[1])
                      for i in range(cnn.image_shape[3][0])]
                     for j in range(nbr_trainer)]
            errC2 = errc2*cnn.dfactive(np.array(filter_weighted_input[1]))
            #error before sigmoid
            errk2 = cnn.kernel_error(np.array(filter_outputs[1]), errC2,
                                     backprop_padding[1])
            if np.amax(np.abs(errk2)) == 0.0:
                print('k2, Gradient Lost')
                Bforbreak = 1
                break
            # 1st convolutional layer
            errF1 = cnn.kernel_backprop_error(errC2, np.array(cnn.filter_weight[1]), 4)
                    # unpool
            errc1 = [[cnn.upsample(errF1[j][i], cnn.image_shape[1][1], pool_structure[0])
                      for i in range(cnn.image_shape[1][0])]
                     for j in range(nbr_trainer)]
            errC1 = errc1*cnn.dfactive(np.array(filter_weighted_input[0]))
            #error before sigmoid
            errk1 = cnn.kernel_error(np.array(filter_outputs[0]), errC1,
                                     backprop_padding[0])
            if np.amax(np.abs(errk1)) == 0.0:
                print('k1, Gradient Lost')
                Bforbreak = 1
                break
            # -------------------------------------------
            # ---------- Weight and Bias shift ----------
            # -------------------------------------------
            
            # output layer
            cnn.neural_weight[1] -= dl * errL.T.dot(al)
            cnn.neural_bias[1] -= dl * np.mean(errL, 0)
            # hidden layer
            cnn.neural_weight[0] -= dl * errl.T.dot(features_vector)
            cnn.neural_bias[0] -= dl * np.mean(errl, 0)
            # convolution 2nd layer
            cnn.filter_weight[1] -= dl * np.array(errk2)
            cnn.filter_bias[1] -= dl * np.einsum('lii', np.mean(errC2, 0))
            # convolution 1st layer
            cnn.filter_weight[0] -= dl * np.array(errk1)
            cnn.filter_bias[0] -= dl * np.einsum('lii', np.mean(errC1, 0))
            
            
            # Prepare for next step
            t += 1
            dl = learning_rate / (1.0 + learning_rate_decay*t)
            b += 1            

        if Bforbreak == 1:
            break
        print('')
        epoch += 1
        if epoch == Nepoch:
            more_epoch = input("More Epoch? (y/n) ")
            if more_epoch == 'y':
                Nepoch += int(input("How many? (int) "))
    step_nbr_record += [t]
    Costs += [Cn]
    
    J +=1
    if J == cnn.nbr_jumps:
        cnn.jump()
    else:
        print('Jump', J)
    
###############################################################################
#                                 RECORDING
###############################################################################
print('')
print('----------------------------RECORDING----------------------------')
index = np.argmin(Costs)*4
sio.savemat('WandB.mat', {'epoch':epoch+epoch_nbr,
                          'filter_weights':cnn.filter_weight,
                          'filter_biases':cnn.filter_bias,
                          'image_shape':cnn.image_shape,
                          'neural_structure':cnn.neural_structure,
                          'neural_weight':cnn.neural_weight,
                          'neural_biases':cnn.neural_bias})
sio.savemat('history_me1.mat', {'acc':Acc, 'cost':Crecord})

##########################################################################
#                                 EVALUATING
##########################################################################
print('')
print('-----------------------------RESULTS-----------------------------')
eval_accuracy = input("Accuracy? (y/n) ")
if eval_accuracy == 'y':
    baches_test = np.arange(nbr_bach_test)
    bb = baches_test*bach_size
    classification = []
    for b in baches_test:
        print('\r' '  bach =', b+1, '/', nbr_bach_test)
        test_data = testerX[bb[b]:bb[b]+bach_size]
        dum, dum, dum, classif, dum, dum = cnn.forward(test_data)
        classification += classif
    acc_test = cnn.accuracy(classification, testerY)
    print(' > accuracy =', "%.2f" % acc_test)

sio.savemat('history_me.mat', {'acc':Acc, 'acc_test':acc_test})
