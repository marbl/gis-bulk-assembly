
#
#  A rather boring test of python.  Might be useful when things break.
#
sub emitPythonTest () {
  print CMD "echo 'PYTHON VERSION:'\n";
  print CMD "command -v python\n";
  print CMD "python --version\n";
  print CMD "\n";
  print CMD "echo ''\n";
  print CMD "echo 'PYTHON LIBRARIES:'\n";
  print CMD "echo  > x.py 'from importlib import metadata'\n";
  print CMD "echo >> x.py 'for dist in metadata.distributions():'\n";
  print CMD "echo >> x.py '  print(dist._path)'\n";
  print CMD "python x.py\n";
  print CMD "\n";
  print CMD "echo ''\n";
  print CMD "echo 'NETWORKX VERSION'\n";
  print CMD "echo  > x.py 'from importlib.metadata import version'\n";
  print CMD "echo >> x.py 'version(\"networkx\")'\n";
  print CMD "python x.py\n";
  print CMD "\n";
  print CMD "echo ''\n";
  print CMD "echo 'NETWORKX TEST'\n";
  print CMD "echo  > x.py 'import networkx as nx'\n";
  print CMD "echo >> x.py 'G = nx.Graph()'\n";
  print CMD "echo >> x.py 'print(G)'\n";
  print CMD "python x.py\n";
  print CMD "\n";
  print CMD "rm x.py\n";
  print CMD "\n";
}


