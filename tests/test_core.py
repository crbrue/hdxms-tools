import pandas as pd
from hdxms.consensus import consensus_from_dataframe
from hdxms.colors import bin_name

def test_consensus_overlap():
    df=pd.DataFrame({'chain':['',''],'start':[1,2],'end':[3,4],'value':[10.0,20.0]})
    out=consensus_from_dataframe(df,mincov=1,value_name='x')
    assert out.loc[out.resi==1,'x'].iloc[0] == 10.0
    assert out.loc[out.resi==2,'x'].iloc[0] == 15.0
    assert out.loc[out.resi==2,'coverage'].iloc[0] == 2

def test_bins():
    assert bin_name(-31)=='hdx_blue3'
    assert bin_name(0)=='hdx_gray0'
    assert bin_name(31)=='hdx_red3'


def test_illustrator_full_and_zoom_exports(tmp_path):
    import pandas as pd
    from hdxms.strip import (
        export_consensus_full_manual_fixed_size,
        export_consensus_zoom_manual_fixed_size,
    )

    csv = tmp_path / "consensus.csv"
    pd.DataFrame({"resi": [1, 2, 3], "pctD_diff": [-20.0, float("nan"), 20.0]}).to_csv(csv, index=False)

    full = tmp_path / "full.svg"
    export_consensus_full_manual_fixed_size(str(csv), str(full))
    text = full.read_text()
    assert 'width="7.4558in"' in text
    assert 'height="0.0751in"' in text
    assert 'preserveAspectRatio="none"' in text
    assert 'fill="#ffffff" stroke="none"' in text
    assert 'fill="#0080FF"' in text
    assert 'fill="#FF0000"' in text

    zoom = tmp_path / "zoom.svg"
    export_consensus_zoom_manual_fixed_size(str(csv), 2, 3, str(zoom))
    text = zoom.read_text()
    assert 'width="1.6089in"' in text
    assert 'height="0.0392in"' in text
    assert 'viewBox="0 0 2 1"' in text
