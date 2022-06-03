set -e # cancels script if any errors occur

# make conda commands accessible in this script
source $(dirname $(which conda))/../etc/profile.d/conda.sh

ENV='splicing_comparison'

# activate env if exists
ENV_INSTALLED=$(conda info --envs | grep "$ENV " || true;)

if [ "$ENV_INSTALLED" ]; then
    conda activate $ENV
fi