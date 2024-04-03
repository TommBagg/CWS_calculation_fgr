# data_preparation_fgr
Algorithms to derive input dataset from LiDAR data and to run forest wind susceptibility calculation in fgr

Author details: Tommaso Baggio, Maximiliano Costa, Niccolò Marchi. Email contact: tommaso.baggio@unipd.it

Copyright statement: This script is the product of the work of Tommaso Baggio, Maximiliano Costa, Niccolò Marchi.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

The release is composed of two algorithms in R environmnet:
- "01_tree_file_V1.R": starting from a CHM raster layer (Canopy Height Model), the script extracts the location of the trees, assigning the height, crown extension, species (using an input shapefile provided by the user). The output is a text file representing the tree location and the relative charcteristics.

- "02_CWS_fgr_DBH_height_crown_gap.R": using the dataset previously created, the script computes the stem density, dominant height, distance to the closest gap and relative extension. Such data, together with the single tree characteristics, are directly passed to fgr for computing the critical wind speeds of overturning and breakage.

The implementation of the algorithms allow the final user to derive the critical wind speeds over large areas (regional scale).

The zip folder provides a training dataset to effectively test the scripts with real data.
