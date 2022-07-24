import numpy as np
import tensorflow as tf
from tensorflow.keras import models, layers

from google.colab import drive
import os
gdrive_path='/content/gdrive/MyDrive/Predicting_interactions_of_m6A-regulated_RNA-binding_proteins/data'
drive.mount('/content/gdrive', force_remount=True)
os.chdir(gdrive_path)
# print(sorted(os.listdir()))

base2int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def load_fasta_methylation_site(fasta):
    with open(fasta) as f:
        for line in f:
            if line[0] != '>':
                raise ValueError(f'Expected FASTA header, got \'{line.strip()[:10]}\'')

            label = int(line.strip()[1])
            
            sequence = f.readline().strip()
            sequence_int = [base2int.get(base, 9999) for base in sequence]
            sequence_onehot = tf.one_hot(sequence_int, depth=4)

            m6A_binding = f.readline().strip()
            m6A_binding_int = tf.constant([int(x) for x in m6A_binding])
            m6A_binding_int = tf.reshape(m6A_binding_int, (1, 200))
            m6A_binding_int = tf.cast( m6A_binding_int, tf.float32)

            numpy_onehot = sequence_onehot.numpy() 
            numpy_m6A = m6A_binding_int.numpy()
            numpy_sequence_m6A = np.concatenate((numpy_onehot, numpy_m6A.T), axis=1)
            sequence_m6A = tf.convert_to_tensor(numpy_sequence_m6A, np.float32)

            yield sequence_m6A, label

def load_fasta_baseline(fasta):
    with open(fasta) as f:
        for line in f:
            if line[0] != '>':
                raise ValueError(f'Expected FASTA header, got \'{line.strip()[:10]}\'')

            label = int(line.strip()[1])

            sequence = f.readline().strip()
            sequence_int = [base2int.get(base, 9999) for base in sequence]
            sequence_onehot = tf.one_hot(sequence_int, depth=4)
    
            yield sequence_onehot, label
            
dataset = tf.data.Dataset.from_generator(lambda: load_fasta_baseline('./IGF2BPs/IGF2BPs_baseline.fasta'), output_signature=(tf.TensorSpec(shape=(200,4), dtype=tf.float32), tf.TensorSpec(shape=(), dtype=tf.float32)))
dataset = dataset.shuffle(len(list(dataset)))
dataset = dataset.cache()         

n_samples = [i for i, _ in enumerate(dataset, start=1)][-1]
dataset_train = dataset.take(int(n_samples * 0.80))
dataset_val = dataset.skip(int(n_samples * 0.80)).take(int(n_samples * 0.1))
dataset_test = dataset.skip(int(n_samples * 0.80) + int(n_samples * 0.1))

dataset_train = dataset_train.batch(128)
dataset_val = dataset_val.batch(128)

def construct_model_baseline():
    model = models.Sequential(name='ML4RG_SS22_Model_BASELINE')

    # input shape (baseline) 150x4
    model.add(layers.Input(shape=(200, 4)))

    model.add(layers.Convolution1D(64, kernel_size=9, activation='relu'))
    model.add(layers.Dropout(0.5, seed = 123))
    model.add(layers.MaxPool1D(2))

    model.add(layers.Convolution1D(64, kernel_size=4, activation='relu'))
    model.add(layers.Dropout(0.5, seed = 123))
    model.add(layers.MaxPool1D(2))

    model.add(layers.Convolution1D(64, kernel_size=4, activation='relu'))
    model.add(layers.Dropout(0.5, seed = 123))
    model.add(layers.MaxPool1D(2))

    model.add(layers.Flatten())
    model.add(layers.Dense(128, activation='relu'))
    model.add(layers.Dropout(0.5, seed = 123))
    model.add(layers.Dense(1, activation='sigmoid'))

    model.compile(loss='binary_crossentropy', optimizer='Adam', metrics=['accuracy'])

    return model
    
model = construct_model() 
history = model.fit(dataset_train, validation_data=dataset_val, epochs = 50)
