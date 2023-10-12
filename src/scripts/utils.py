import pandas as pd
import re


def replace_colnames(df):

    # Extract only those between ERR and _
    for col in df.columns:
        if "ERR" in col:
            new_col = re.search('ERR(.*)_', col).group(1)
            df.rename(columns={col: str(new_col)}, inplace=True)

    return df


def txt_2_dataframe(txt_filePath):

    df = pd.read_csv(f'{txt_filePath}', delimiter="\t")
    return df

# data_renamed["Sum"] = data_renamed.sum(axis=1)
# data_renamed = data_renamed.sort_values(by=["Sum"], ascending=False)
