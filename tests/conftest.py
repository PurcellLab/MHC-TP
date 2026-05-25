import textwrap
import pytest


@pytest.fixture
def pssm_text(tmp_path):
    """A 2-position PSSM file in the NetMHCpan/Gibbs matrix layout."""
    content = textwrap.dedent("""\
        # comment line
        Pos A R N D C Q E G H I L K M F P S T W Y V
        1 P 0.1 0.2 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.7
        2 P 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
        """)
    p = tmp_path / "mat.txt"
    p.write_text(content)
    return p
