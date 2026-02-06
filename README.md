<h1>ORF Gene Reader</h1>
<h3>About</h3>
<text>This is a python based bioinformatics tool for analysing DNA sequences from FASTA files and identifying open reading frames (ORFs) across all six possible reading frames.</text>
<h3>Goal</h3>
<ol>
  <li>Read a DNA sequence from FASTA file</li>
  <li>Generate the complementary and reverse complementary strands</li>
  <li>Translate codons into amino acids</li>
  <li>Display amino acid sequences for all six reading frames</li>
</ol>
<h3>How it Works</h3>
<h4>1. FASTA Input</h4>

<p>Run the script using:</p>

<pre><code class="language-bash">python readORF.py</code></pre>

<p>Then select the FASTA file</p>
<img src="ORFGeneReader/fasta.PNG" alt="fasta file" width="300">
<p>The program reads the a FASTA file and extracts</p>
<ul>
  <li>Description/header line</li>
  <li>DNA sequence data</li>
  <li>Newline characters(for continous nucleotide sequence)</li>
</ul>

<h4>2. Complementary and Reverse Complenetary Strands</h4>
<p>The complementary strand is generated using stabdard base pairing rules</p>
<table>
  <thead>
    <tr>
      <th>Base</th>
      <th>Complement</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>A</td>
      <td>T</td>
    </tr>
    <tr>
      <td>T</td>
      <td>A</td>
    </tr>
    <tr>
      <td>C</td>
      <td>G</td>
    </tr>
    <tr>
      <td>G</td>
      <td>C</td>
    </tr>
  </tbody>
</table>



