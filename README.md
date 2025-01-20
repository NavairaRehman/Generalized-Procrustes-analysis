
# Generalized Procrustes analysis (GPA) for face alignment 
GPA is a method of statistical analysis that can be used to compare the shapes of objects. This piece of code employs GPA for face alignment : developed by Pulak Purkait. The details of the algorithms can be found in https://graphics.stanford.edu/courses/cs164-09-spring/Handouts/paper_shape_spaces_imm403.pdf

For any queries contact:  Email - pulak.isi@gmail.com 

## Usage 

The dataset can be found at https://ibug.doc.ic.ac.uk/resources/facial-point-annotations/

Extract the datasets:

    cat 300w.zip.* > 300w.zip; unzip 300w.zip 

To execute the code:

    run GPA_PCA.m in MATLAB


## 3D Generalized Procrustes Analysis (GPA)

1. **Input Files**: Place your `.obj` files (containing vertices and faces) in the `input` directory.
2. **Run the Script**: Execute the `GPA_3d.m` script in MATLAB. This performs GPA on the vertices while retaining the original faces.
   ```matlab
   GPA_3d
   ```
3. **Output**: Transformed `.obj` files are saved in the `output` directory with the same structure.

Ensure the input files follow the standard OBJ format for proper processing.

