import React, { Component } from "react";
import PropTypes from "prop-types";
import getCookie from "./utils";

function HelpText(props) {
  const iconStyleOverride = {
    float: 'right',
    fontWeight: 'bold',
  };

  return (
    <button type="button" className="button is-small" style={iconStyleOverride} onClick={() => alert(props.text)}>?</button>
  );
}

class JobForm extends Component {
  static propTypes = {
    endpoint: PropTypes.string.isRequired
  };

  constructor(props) {
    super(props);
    this.state = {
      title: "",
      description: "",
      email: "",
      foreground_data: "",
      foregroundformat: "",
      foreground_filename: "",
      contextformat: "2",
      context_data: "",
      context_filename: "",
      p_value_cutoff: 0.001,
      position_specific: true,
      minimum_occurrences: 20,
      fold_change_cutoff: 1.0,
      max_depth: "",
      width: 15,
      center_sequences: false,
      multiple_testing_correction: true,
      positional_weighting: true,
      compound_residues: true,
      compound_residue_decomposition: true,
      compoundresidue_file: "",
      compoundresidue_filename: "",
      require_context_id: true,
      extension_direction: "1",
      redundancylevel: "1",
      first_protein_only: true,
      originalrowmerge: "3",
      cluster_sequences: true,
    };
    this.toggleForegroundFormat = this.toggleForegroundFormat.bind(this);
    this.handleChange = this.handleChange.bind(this);
    this.handleCheckboxChange = this.handleCheckboxChange.bind(this);
    this.toggleAdvancedOptions = this.toggleAdvancedOptions.bind(this);
    this.handleUpload = this.handleUpload.bind(this);
    this.handleSubmit = this.handleSubmit.bind(this);
    this.resetForm = this.resetForm.bind(this);
  }

  toggleForegroundFormat = () => {
    this.setState({

    });
  }

  handleChange = e => {
    this.setState({ [e.target.name]: e.target.value });
    if (e.target.name == "foregroundformat") {
      if (e.target.value == 1) {
        this.setState({ redundancylevel: 1 });
      }
      else if (e.target.value == 3) {
        this.setState({ redundancylevel: 2 });
      }
    }
    else if (e.target.name == "width") {
      this.setState({ center_sequences: (e.target.value % 2 == 0) ? true : false });
    }
  };

  handleCheckboxChange = e => {
    this.setState({ [e.target.name]: e.target.checked });
  };

  toggleAdvancedOptions = () => {
    var advancedOptions = document.getElementById('advancedOptions');
    if (advancedOptions.style.display === "none") {
      document.getElementById('advanced-options-header').innerHTML = 'Advanced options &#9660';
      advancedOptions.style.display = "block";
    } else {
      document.getElementById('advanced-options-header').innerHTML = 'Advanced options &#9650';
      advancedOptions.style.display = "none";
    }
  }

  handleUpload = e => {
    this.setState({ [e.target.name]: e.target.files[0] });
    this.setState({
      [e.target.name.slice(0, e.target.name.indexOf("_")) + "_filename"]: e.target.files[0].name
    });
    document.getElementById(e.target.name + "-text").textContent = e.target.files[0].name;
  };
  
  resetForm = () => {
    {/*
    if (document.getElementById('peptidelist').checked) {
      document.getElementById('peptidelist').checked = false;
    }
    if (document.getElementById('prealigned').checked) {
      document.getElementById('prealigned').checked = false;
    }
    */}
    
    document.getElementById('advanced-options-header').innerHTML = 'Advanced options &#9650';
    advancedOptions.style.display = "none";

    if (this.state.foreground_data != "") {
      document.getElementById("foreground_data-text").textContent = 'Choose a file...';
    }
    if (this.state.context_data != "") {
      document.getElementById("context_data-text").textContent = 'Choose a file...';
    }

    this.setState({
      title: "",
      email: "",
      description: "",
      foreground_data: "",
      foregroundformat: "",
      foreground_filename: "",
      contextformat: "2",
      context_data: "",
      context_filename: "",
      p_value_cutoff: 0.001,
      contextformat: "",
      position_specific: true,
      minimum_occurrences: 20,
      fold_change_cutoff: 1.0,
      max_depth: "",
      width: 15,
      center_sequences: false,
      center_sequences: false,
      multiple_testing_correction: true,
      positional_weighting: true,
      compound_residues: true,
      compound_residue_decomposition: true,
      compoundresidue_file: "",
      compoundresidue_filename: "",
      require_context_id: true,
      extension_direction: "1",
      redundancylevel: "1",
      first_protein_only: true,
      originalrowmerge: "1",
      cluster_sequences: true,
    });
  }

  handleSubmit = e => {
    e.preventDefault();
    let formData = new FormData();
    const {
      title,
      description,
      email,
      foreground_data,
      foregroundformat,
      foreground_filename,
      contextformat,
      context_data,
      context_filename,
      p_value_cutoff,
      position_specific,
      minimum_occurrences,
      fold_change_cutoff,
      width,
      center_sequences,
      multiple_testing_correction,
      positional_weighting,
      compound_residues,
      compound_residue_decomposition,
      compoundresidue_file,
      compoundresidue_filename,
      require_context_id,
      extension_direction,
      redundancylevel,
      first_protein_only,
      originalrowmerge,
      cluster_sequences
    } = this.state;
    formData.append('foreground_data', foreground_data);
    formData.append('foregroundformat', foregroundformat);
    formData.append('foreground_filename', foreground_filename);
    formData.append('contextformat', contextformat);
    formData.append('context_data', context_data);
    formData.append('context_filename', context_filename);
    formData.append('title', title);
    formData.append('description', description);
    formData.append('email', email);
    formData.append('p_value_cutoff', p_value_cutoff);
    formData.append('position_specific', position_specific);
    formData.append('minimum_occurrences', minimum_occurrences);
    formData.append('fold_change_cutoff', fold_change_cutoff);
    formData.append('width', width);
    formData.append('center_sequences', center_sequences);
    formData.append('multiple_testing_correction', multiple_testing_correction);
    formData.append('positional_weighting', positional_weighting);
    formData.append('compound_residues', compound_residues);
    formData.append('compound_residue_decomposition', compound_residue_decomposition);
    formData.append('compoundresidue_file', compoundresidue_file);
    formData.append('compoundresidue_filename', compoundresidue_filename);
    formData.append('require_context_id', require_context_id);
    formData.append('extension_direction', extension_direction);
    formData.append('redundancylevel', redundancylevel);
    formData.append('first_protein_only', first_protein_only);
    formData.append('originalrowmerge', originalrowmerge);
    formData.append('cluster_sequences', cluster_sequences);

    const csrftoken = getCookie('csrftoken');
    const conf = {
      method: "post",
      body: formData,
      headers: new Headers({ "X-CSRFToken": csrftoken })
    };

    fetch(this.props.endpoint, conf).then(response => console.log(response));
    
    {/* this.resetForm(); */}

    alert("Thank you for your submission. Your results will be emailed to you.");
  };
  
  render() {
    const {
      title,
      description,
      email,
      foreground_data,
      foregroundformat,
      foreground_filename,
      contextformat,
      context_data,
      context_filename,
      p_value_cutoff,
      position_specific,
      minimum_occurrences,
      fold_change_cutoff,
      width,
      center_sequences,
      multiple_testing_correction,
      positional_weighting,
      compound_residues,
      compound_residue_decomposition,
      compoundresidue_file,
      compoundresidue_filename,
      require_context_id,
      extension_direction,
      redundancylevel,
      first_protein_only,
      originalrowmerge,
      cluster_sequences
    } = this.state;
    return (
      <section className="jobformcontainer">
        <h4 className="title is-4">Submit a new job for analysis.</h4>
        <form className="jobform" onSubmit={this.handleSubmit} encType="multipart/form-data">
          <div className="field">
            <label className="label">Email</label>
            <div className="control">
              <input
                className="input"
                type="email"
                name="email"
                onChange={this.handleChange}
                value={email}
                required
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Title (20 characters max.)</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="title"
                maxlength="20"
                onChange={this.handleChange}
                value={title}
                required
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Description</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="description"
                onChange={this.handleChange}
                value={description}
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Select foreground data file to upload.</label>
            <label className="file-label">
              <div className="control">
                <input
                  className="file-input"
                  type="file"
                  name="foreground_data"
                  onChange={this.handleUpload}
                  autoComplete="off"
                  required
                />
                <span className="file-cta">
                  <span className="file-icon">
                    <i className="fas fa-upload"></i>
                  </span>
                  <span className="file-label" id="foreground_data-text">
                    Choose a file…
                  </span>
                </span>
              </div>
            </label>
            <HelpText text={
                  "Column names may not contain special characters other than underscores.\n\n"
                  + "Cell values of T and F in extra data columns will be interpreted as True and False respectively.\n"} />
          </div>
          <div>
            <label className="label">Select foreground data set format.</label>
            <div>
              <label>
                <input
                  type="radio"
                  name="foregroundformat"
                  id="prealigned"
                  onChange={this.handleChange}
                  value="1"
                  checked={foregroundformat == "1"}
                  required
                />
                Prealigned text file (<a href='/rolim/textfile' download>Example</a>)
              </label>
            </div>
            <div>
              <label>
                <input
                  type="radio"
                  name="foregroundformat"
                  id="peptidelist"
                  onChange={this.handleChange}
                  value="3"
                  checked={foregroundformat == "3"}
                />
                Text file peptide list (<a href='/rolim/peptidelist' download>Example</a>)
              </label>
            </div>
          </div>
          <br />
          <div>
            <label className="label">Select context data set format.</label>
            <div>
              <label>
                <input
                  type="radio"
                  name="contextformat"
                  id="swissprot-human"
                  onChange={this.handleChange}
                  value="2"
                  checked={contextformat == "2"}
                  required
                />
                Swiss-Prot Human
              </label>
            </div>
            {/*
            <div>
              <label>
                <input
                  type="radio"
                  name="contextformat"
                  id="swissprot-mouse"
                  onChange={this.handleChange}
                  value="3"
                  checked={contextformat == "3"}
                />
                Swiss-Prot Mouse
              </label>
            </div>
            */}
            <div>
              <label>
                <input
                  type="radio"
                  name="contextformat"
                  id="fasta"
                  onChange={this.handleChange}
                  value="1"
                  checked={contextformat == 1}
                />
                Other (uploaded FASTA file)
              </label>
            </div>
          </div>
          <br />
          <div className="field">
            <label className="label">Select context FASTA file to upload. (Optional)</label>
            <label className="file-label">
              <div className="control">
                <input
                  className="file-input"
                  type="file"
                  name="context_data"
                  onChange={this.handleUpload}
                  autoComplete="off"
                />
                <span className="file-cta">
                  <span className="file-icon">
                    <i className="fas fa-upload"></i>
                  </span>
                  <span className="file-label" id="context_data-text">
                    Choose a file…
                  </span>
                </span>
              </div>
            </label>
          </div>
          <br />
          <div className="field">
              <label className="checkbox">
                Enter desired width of expanded sequences. (MEROPS comparison only supported for width-8 sequences)
              </label>
              <HelpText text={
                  "The number of residues in each aligned sequence of the foreground data set. If peptides are"
                  + " supplied and extension/alignment is enabled, this is the final length of each extended and aligned sequence."
                  + " If pre-aligned sequences are supplied, each supplied sequence MUST be of this length.\n\n"
                  + "*Width must be set to 8 to enable enhanced protease analysis.\n"} />
              <div className="control">
                <input
                  className="input"
                  type="number"
                  name="width"
                  onChange={this.handleChange}
                  value={width}
                  min="1"
                  max="25"
                />
              </div>
          </div>
          <br />
          <h5 id="advanced-options-header" className="title is-5" onClick={this.toggleAdvancedOptions}>Advanced options &#9650;</h5>
          <br />
          <div id="advancedOptions" style={{display: "none"}}>
            <div className="field">
              <div className="control">
                <label className="checkbox">
                  <input
                    type="checkbox"
                    name="center_sequences"
                    onChange={this.handleCheckboxChange}
                    value={center_sequences}
                    checked={center_sequences}
                  />
                  Center sequence position numbers?
                </label>
                <HelpText text={
                    "For example:\n\nCentered:\n\n"
                    + "    p4-p3-p2-p1-p1'-p2'-p3'-p4\n\n"
                    + "Non-centered:\n\n"
                    + "    p1-p2-p3-p4-p5-p6-p7-p8-p9-p10-p11-p12-p13-p14-p15\n\n"
                    + "*Protease analysis only supported for centered sequences with a window width of 8.\n"} />
              </div>
            </div>
            <div className="field">
              <div className="control">
                <label className="checkbox">
                  <input
                    type="checkbox"
                    name="compound_residues"
                    onChange={this.handleCheckboxChange}
                    value={compound_residues}
                    checked={compound_residues}
                  />
                Detect compound residue groups?
                </label>
                <HelpText text={
                    "Enables aggregation of single amino"
                    + " acids into groups of biochemically and/or structurally similar amino acids which"
                    + " may be cumulatively enriched.\n"} />
              </div>
            </div>
            <div className="field">
              <div className="control">
                <label className="checkbox">
                  <input
                    type="checkbox"
                    name="compound_residue_decomposition"
                    onChange={this.handleCheckboxChange}
                    value={compound_residue_decomposition}
                    checked={compound_residue_decomposition}
                  />
                  Enable compound residue decomposition?
                </label>
                <HelpText text={
                    "Enables decomposition of enriched"
                    + " compound positional residue groups into subsets composed of the constituents of the compound"
                    + " residue group (e.g. [RK] -> [R, K])\n"} />
              </div>
            </div>
            <div className="field">
            <label className="label">Select custom compound residue file to upload.</label>
            <label className="file-label">
              <div className="control">
                <input
                  className="file-input"
                  type="file"
                  name="compoundresidue_file"
                  onChange={this.handleUpload}
                  autoComplete="off"
                />
                <span className="file-cta">
                  <span className="file-icon">
                    <i className="fas fa-upload"></i>
                  </span>
                  <span className="file-label" id="compoundresidue_file-text">
                    Choose a file…
                  </span>
                </span>
                <a href='/rolim/compoundresiduefile' download>Example</a>
              </div>
            </label>
          </div>
            <div className="field">
              <div className="control">
                <label className="checkbox">  
                  <input
                    type="checkbox"
                    name="multiple_testing_correction"
                    onChange={this.handleCheckboxChange}
                    value={multiple_testing_correction}
                    checked={multiple_testing_correction}
                  />
                  Enable multiple testing correction?
                </label>
                <HelpText text={
                    "Enables optional Bonferroni correction"
                    + " for positional residue p-values.\n"} />
              </div>
            </div>
            <div className="field">
              <div className="control">
                <label className="checkbox">
                  <input
                    type="checkbox"
                    name="positional_weighting"
                    onChange={this.handleCheckboxChange}
                    value={positional_weighting}
                    checked={positional_weighting}
                  />
                  Enable positional weighting?
                </label>
                <HelpText text={
                    "Enables optional positional weighting term in"
                    + " positional residue enrichment calculation. Positional weight is calculated as"
                    + " (1 / # distinct residues in a position).\n"} />
              </div>
            </div>
            <div className="field">
              <div className="control">
                <label className="checkbox">
                  <input
                    type="checkbox"
                    name="position_specific"
                    onChange={this.handleCheckboxChange}
                    value={position_specific}
                    checked={position_specific}
                  />
                  Position specific background?
                </label>
                <HelpText text={
                    "Enables background frequency calculation"
                    + " from a complete, position-specific background derived from the context data set used for"
                    + " an analysis. When disabled, background frequencies are averaged across all positions of the"
                    + " context data set and dynamically updated when position/residue pairs are eliminated from"
                    + " the foreground data set.\n"} />
              </div>
            </div>
            <div className="field">
              <div className="control">
                <label className="checkbox">
                  <input
                    type="checkbox"
                    name="require_context_id"
                    onChange={this.handleCheckboxChange}
                    value={require_context_id}
                    checked={require_context_id}
                  />
                  Require protein identifier?
                </label>
                <HelpText text={
                    "Require a protein identifier for each foreground"
                    + " sequence. Foreground protein identifiers must match the format of protein identifiers used"
                    + " in the context data set.\n"} />
              </div>
            </div>
            <div>
              <label className="label">Select sequence extension direction.</label>
              <HelpText text={
                  "This tool supports alignment and extension of unaligned foreground data sets using the selected context data set.\n\n"
                  + "N-terminal: Supplied sequences will be extended in the N-terminal direction up to the specified aligned sequence width (default).\n"
                  + "C-terminal: Supplied sequences will be extended in the C-terminal direction up to the specified aligned sequence width.\n"
                  + "Both: Supplied sequences will be extended in both the N-terminal and C-terminal directions up to the specified aligned sequence width,"
                  + " producing a separate sequence for each extension direction.\n"} />
              <div>
                <label>
                  <input
                    type="radio"
                    name="extension_direction"
                    id="nextension"
                    onChange={this.handleChange}
                    value="1"
                    checked={extension_direction == 1}
                    required
                  />
                  N-terminal
                </label>
              </div>
              <div>
                <label>
                  <input
                    type="radio"
                    name="extension_direction"
                    id="cextension"
                    onChange={this.handleChange}
                    value="2"
                    checked={extension_direction == 2}
                  />
                  C-terminal
                </label>
              </div>
              <div>
                <label>
                  <input
                    type="radio"
                    name="extension_direction"
                    id="bothextension"
                    onChange={this.handleChange}
                    value="3"
                    checked={extension_direction == 3}
                  />
                  Both
                </label>
            </div>
          </div>
          <br />
          <div className="field">
            <label className="checkbox">
              P-value threshold.
            </label>
            <HelpText text={
                "The p-value corresponding to the frequency of a position/residue pair must be below this threshold"
                + " in order to be considered significantly enriched. This is determined based on the iteration-specific"
                + " foreground/background context in which the pair is observed.\n\n"
                + "Note that this is not the same as the p-value corresponding to the frequency of a motif in your foreground"
                + " data set. RoLiM does not calculate the p-values of motif foreground frequencies, and the final motif list"
                + " returned by RoLiM is not filtered based on motif-level significance.\n\n"} />
            <div className="control">
              <input
                className="input"
                type="number"
                name="p_value_cutoff"
                onChange={this.handleChange}
                value={p_value_cutoff}
                min="0"
                max="1"
                step="any"
              />
            </div>
          </div>
          <div className="field">
            <label className="checkbox">
              Minimum occurrences
            </label>
            <HelpText text={
                "The minimum frequency of a position/residue pair in the foreground data set"
                + " required for the pair to be considered enriched.\n"} />
            <div className="control">
              <input
                className="input"
                type="number"
                name="minimum_occurrences"
                onChange={this.handleChange}
                value={minimum_occurrences}
                min="1"
              />
            </div>
          </div>
            {/* <div className="field">
              <label className="checkbox">
                Fold difference cutoff
              </label>
              <HelpText text={
                  "The minimum fold difference of position/residue pair in the foreground data set"
                  + " vs. the background data set required for the pair to be considered enriched.\n"} />
              <div className="control">
                <input
                  className="input"
                  type="number"
                  name="fold_change_cutoff"
                  onChange={this.handleChange}
                  value={fold_change_cutoff}
                  min="1"
                  step="any"
                />
              </div>
            </div> */}
            <div className="field">
              <div className="control">
                <label className="checkbox">
                  <input
                    type="checkbox"
                    name="first_protein_only"
                    onChange={this.handleCheckboxChange}
                    value={first_protein_only}
                    checked={first_protein_only}
                  />
                  Only use first protein identifier provided for each peptide?
                </label>
                <HelpText text={
                    "If selected, only the first protein identifer for each row in the foreground"
                    + " data set will be used for alignment and extenson. Otherwise, all provided"
                    + " protein identifiers for each row will be used.\n"} />
              </div>
            </div>
            <div>
                <label className="label">Select sequence redundancy elimination level.</label>
                <HelpText text={
                    "This tool supports mulitple levels of sequence redundancy elimination.\n\n"
                    + "None: No sequence redundancy elimination will be performed (default for pre-aligned foreground data sets).\n"
                    + "Protein: Redundant sequences from the same protein position will be eliminated (default for unaligned foreground data sets).\n"
                    + "Sequence: All redundant sequences will be eliminated.\n"} />
                <div>
                  <label>
                    <input
                      type="radio"
                      name="redundancylevel"
                      id="none"
                      onChange={this.handleChange}
                      value="1"
                      checked={redundancylevel == 1}
                      required
                    />
                    None
                  </label>
                </div>
                <div>
                  <label>
                    <input
                      type="radio"
                      name="redundancylevel"
                      id="protein"
                      onChange={this.handleChange}
                      value="2"
                      checked={redundancylevel == 2}
                    />
                    Protein
                  </label>
                </div>
                <div>
                  <label>
                    <input
                      type="radio"
                      name="redundancylevel"
                      id="sequence"
                      onChange={this.handleChange}
                      value="3"
                      checked={redundancylevel == 3}
                    />
                    Sequence
                  </label>
              </div>
            </div>
            <br />
            <div>
                <label className="label">Select level at which to merge multiple peptide alignments.</label>
                <HelpText text={
                    "Peptides may match more than one position in the context proteome. This setting specifies how to handle those multiple matches.\n\n"
                    + "None: Multiple matches will not be merged. Each unique match will be included as a separate sequence in the aligned foreground data set.\n"
                    + "Protein: Multiple matches mapping to the same protein will be merged, replacing positions in the extended region of the sequence which disagree with an X.\n"
                    + "All: All multiple matches will be merged resulting in one aligned sequence, possibly containing Xs in the extended region of the aligned sequence (default for unaligned foreground data sets).\n"} />
                <div>
                  <label>
                    <input
                      type="radio"
                      name="originalrowmerge"
                      id="none"
                      onChange={this.handleChange}
                      value="1"
                      checked={originalrowmerge == 1}
                      required
                    />
                    None
                  </label>
                </div>
                <div>
                  <label>
                    <input
                      type="radio"
                      name="originalrowmerge"
                      id="protein"
                      onChange={this.handleChange}
                      value="2"
                      checked={originalrowmerge == 2}
                    />
                    Protein
                  </label>
                </div>
                <div>
                  <label>
                    <input
                      type="radio"
                      name="originalrowmerge"
                      id="all"
                      onChange={this.handleChange}
                      value="3"
                      checked={originalrowmerge == 3}
                    />
                    All
                  </label>
              </div>
            </div>
            <br />
            <div className="field">
              <div className="control">
                <label className="checkbox">
                  <input
                    type="checkbox"
                    name="cluster_sequences"
                    onChange={this.handleCheckboxChange}
                    value={cluster_sequences}
                    checked={cluster_sequences}
                  />
                  Enable hierarchical clustering of sequences and motifs?
                </label>
                <HelpText text={
                    "Enables optional hierarchical clustering of sequences and motifs.\n\n"
                    + " This is a time- and memory-intensive operation, and should\n"
                    + " be disabled for very large foreground data sets.\n"} />
              </div>
            </div>
          </div>
          <br />
          <div className="control">
            <button type="submit" className="button is-info">
              Submit job
            </button>
          </div>
        </form>
      </section>
    );
  }
}

export default JobForm;