library("tmod")
# A Small Example (Boston Housing Data)

# Building a model in Keras starts by constructing an empty Sequential model.

load("~/Koop_Domaszewska/Data/MDS.RDa")
load("~/Koop_Domaszewska/Data/data_matrix.RDa")
load("~/Koop_Domaszewska/Data/myIFN_I_set.RDa")
load("~/Koop_Domaszewska/Data/study_donor_vector.RDa")

#meta_info = rf_meta$targets
#write.table(meta_info,"~/Koop_Domaszewska/Misc/Meta_Information.tsv",sep ="\t", row.names =T, quote = F)
meta_info = read.table("~/Koop_Domaszewska/Misc/Meta_Information.tsv", sep ="\t", header = T, stringsAsFactors = F, row.names = 1)
meta_match = match(colnames(data_matrix), rownames(meta_info), nomatch = 0)
meta_data = meta_info[meta_match,]
dim(data_matrix)
data_matrix = data_matrix[ ,meta_match != 0]
dim(data_matrix)

data_train = data_matrix[, meta_data$study != "bloom"]
data_test = data_matrix[, meta_data$study == "bloom"]

label_train = meta_data$activeTB[ meta_data$study != "bloom"]
label_test  = meta_data$activeTB[ meta_data$study == "bloom"]

dim(data_train)
dim(data_test)
length(label_train)
length(label_test)

#Sys.setenv(OPENBLAS_NUM_THREADS="16")
#Sys.setenv(NUMEXPR_NUM_THREADS="16")
#Sys.setenv(OMP_NUM_THREADS="16")

# below commands will throw error don't panic!
library(tensorflow) # only loaded earlier for setting parameters
#config = tf$ConfigProto()
#result = tryCatch({ config$gpu_options$allow_growth=TRUE }, error = function(e) { })
#print(config$gpu_options$allow_growth)
#config$inter_op_parallelism_threads = 1L
#config$intra_op_parallelism_threads = 1L
#server = tf$train$Server$create_local_server(config=config)
#sess = tf$Session(server$target)

# load libraries
library(kerasR) # if you get an error message -> install.packages("kerasR")
library(keras)

#k_set_session(sess)

nr_epochs = 20 # too low in practice but needed due to performance issues

mod <- Sequential()
# The result of Sequential, as with most of the functions provided by kerasR, is a python.builtin.object. This object type, defined from the reticulate package, provides direct access to all of the methods and attributes exposed by the underlying python class. To access these, we use the $ operator followed by the method name. Layers are added by calling the method add. This function takes as an input another python.builtin.object, generally constructed as the output of another kerasR function. For example, to add a dense layer to our model we do the following:

mod$add(Dense(units = 50, input_shape = nrow(data_train)))
# We have now added a dense layer with 200 neurons. The first layer must include a specification of the input_shape, giving the dimensionality of the input data. Here we set the number of input variables equal to 13. Next in the model, we add an activation defined by a rectified linear unit to the model:

mod$add(Activation("relu"))
mod$add(Dense(units = 1))

# Once the model is fully defined, we have to compile it before fitting its parameters or using it for prediction. Compiling a model can be done with the method compile, but some optional arguments to it can cause trouble when converting from R types so we provide a custom wrapper keras_compile. At a minimum we need to specify the loss function and the optimizer. The loss can be specified with just a string, but we will pass the output of another kerasR function as the optimizer. Here we use the RMSprop optimizer as it generally gives fairly good performance:

keras_compile(mod,  loss = 'mse', optimizer = RMSprop())
# Now we are able to fit the weights in the model from some training data, but we do not yet have any data from which to train! Let’s load some using the wrapper function load_boston_housing. We provide several data loading functions as part of the package, and all return data in the same format. In this case it will be helpful to scale the data matrices:

#boston <- load_boston_housing()
X_train = scale( data_train )
Y_train <- label_train
Y_train[Y_train == "TB"] = 1.0
Y_train[Y_train != 1.0 ] = 0.0
Y_train = as.integer(Y_train)
X_test <- scale(  data_test )
Y_test <- label_test
Y_test[Y_test == "TB"] = 1.0
Y_test[Y_test != 1.0 ] = 0.0
Y_test = as.intger(Y_test)

# Now, we call the wrapper keras_fit in order to fit the model from this data. As with the compilation, there is a direct method for doing this but you will likely run into data type conversion problems calling it directly. Instead, we see how easy it is to use the wrapper function (if you run this yourself, you will see that Keras provides very good verbose output for tracking the fitting of models):

keras_fit(
    mod,
    t(X_train),
    Y_train,
    batch_size = 32,
    epochs = nr_epochs,
    verbose = 1,
    validation_split = 0.1
)

## output

pred = keras_predict(mod, normalize(t(X_test)))
sd(as.numeric(round(pred,0)) - Y_test) / sd(Y_test)

table(cbind(round(pred,1),Y_test))

## A Larger Example (MNIST)
# To show the power of neural networks we need a larger dataset to make use of. A popular first dataset for applying neural networks is the MNIST Handwriting dataset, consisting of small black and white scans of handwritten numeric digits (0-9). The task is to build a classifier that correctly identifies the numeric value from the scan. We may load this dataset in with the following:

#mnist <- load_mnist()
#X_train <- mnist$X_train
#Y_train <- mnist$Y_train
#X_test <- mnist$X_test
#Y_test <- mnist$Y_test
#dim(X_train)

X_train <- array(X_train, dim = c(dim(X_train)[1], prod(dim(X_train)[-1]))) / 255
X_test <- array(X_test, dim = c(dim(X_test)[1], prod(dim(X_test)[-1]))) / 255
# Finally, we want to process the response vector y into a different format as well. By default it is encoded in a one-column matrix with each row giving the number represented by the hand written image. We instead would like this to be converted into a 10-column binary matrix, with exactly one 1 in each row indicating which digit is represented. This is similar to the factor contrasts matrix one would construct when using factors in a linear model. In the neural network literature it is call the one-hot representation. We construct it here via the wrapper function to_categorical. Note that we only want to convert the training data to this format; the test data should remain in its original one-column shape.

Y_train <- to_categorical(mnist$Y_train, 10)
# With the data in hand, we are now ready to construct a neural network. We will create three blocks of identical Dense layers, all having 512 nodes, a leaky rectified linear unit, and drop out. These will be followed on the top output layer of 10 nodes and a final softmax activation. These are fairly well-known choices for a simple dense neural network and allow us to show off many of the possibilities within the kerasR interface:

mod <- Sequential()

mod$add(Dense(units = 512, input_shape = dim(X_train)[1]))
mod$add(LeakyReLU())
mod$add(Dropout(0.25))

mod$add(Dense(units = 512))
mod$add(LeakyReLU())
mod$add(Dropout(0.25))

mod$add(Dense(units = 512))
mod$add(LeakyReLU())
mod$add(Dropout(0.25))

mod$add(Dense(1))
mod$add(Activation("softmax"))
# We then compile the model with the “categorical_crossentropy” loss and fit it on the training data:
# categorical_crossentropy
# sparse_categorical_crossentropy
#keras_compile(mod,  loss = 'categorical_crossentropy', optimizer = RMSprop())
keras_compile(mod,  loss = 'sparse_categorical_crossentropy', optimizer = RMSprop())

Y_train[Y_train==1] = 0.9

keras_fit(
  mod,
  t(X_train),
  Y_train,
  batch_size = 32,
  epochs = nr_epochs,
  verbose = 1,
  validation_split = 0.1
)
# Now that the model is trained, we could use the function keras_predict once again, however this would give us an output matrix with 10 columns. It is not too much work to turn this into predicted classes, but kerasR provides keras_predict_classes that extracts the predicted classes directly. Using this we are able to evaluate the data on the test set.

Y_test_hat <- keras_predict_classes(mod, t(X_test))
table( cbind(Y_test, Y_test_hat))
mean(Y_test == Y_test_hat)

