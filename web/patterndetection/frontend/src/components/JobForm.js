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

  state = {
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
    require_context_id: true,
    extension_direction: "1",
    redundancylevel: "1",
    first_protein_only: true,
    originalrowmerge: "3"
  };

  handleChange = e => {
    this.setState({ [e.target.name]: e.target.value });
    if (e.target.name == "foregroundformat") {
      if (e.target.value == 1) {
        this.setState({ "redundancylevel": 1 });
      }
      else if (e.target.value == 3) {
        this.setState({ "redundancylevel": 2 });
      }
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
  
  handleSubmit = e => {
    e.preventDefault();
    let formData = new FormData();
    formData.append('foreground_data', this.state.foreground_data);
    formData.append('foregroundformat', this.state.foregroundformat);
    formData.append('foreground_filename', this.state.foreground_filename);
    formData.append('contextformat', this.state.contextformat);
    formData.append('context_data', this.state.context_data);
    formData.append('context_filename', this.state.context_filename);
    formData.append('title', this.state.title);
    formData.append('description', this.state.description);
    formData.append('email', this.state.email);
    formData.append('p_value_cutoff', this.state.p_value_cutoff);
    formData.append('position_specific', this.state.position_specific);
    formData.append('minimum_occurrences', this.state.minimum_occurrences);
    formData.append('fold_change_cutoff', this.state.fold_change_cutoff);
    formData.append('width', this.state.width);
    formData.append('center_sequences', this.state.center_sequences);
    formData.append('multiple_testing_correction', this.state.multiple_testing_correction);
    formData.append('positional_weighting', this.state.positional_weighting);
    formData.append('compound_residues', this.state.compound_residues);
    formData.append('compound_residue_decomposition', this.state.compound_residue_decomposition);
    formData.append('require_context_id', this.state.require_context_id);
    formData.append('extension_direction', this.state.extension_direction);
    formData.append('redundancylevel', this.state.redundancylevel);
    formData.append('first_protein_only', this.state.first_protein_only);
    formData.append('originalrowmerge', this.state.originalrowmerge);

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
      require_context_id,
      extension_direction,
      redundancylevel,
      first_protein_only,
      originalrowmerge
    } = this.state;
    const job = {
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
      require_context_id,
      extension_direction,
      redundancylevel,
      first_protein_only,
      originalrowmerge
    };
    const csrftoken = getCookie('csrftoken');
    const conf = {
      method: "post",
      body: formData,
      headers: new Headers({ "X-CSRFToken": csrftoken })
    };

    fetch(this.props.endpoint, conf).then(response => console.log(response));
    
    if (document.getElementById('peptidelist').checked) {
      document.getElementById('peptidelist').checked = false;
    }
    if (document.getElementById('prealigned').checked) {
      document.getElementById('prealigned').checked = false;
    }
    if (this.state.foreground_data != "") {
      document.getElementById("foreground_data-text").textContent = 'Choose a file...';
    }
    if (this.state.foreground_data != "") {
      document.getElementById("context_data-text").textContent = 'Choose a file...';
    }
    this.setState({ ['title']: "" });
    this.setState({ ['email']: "" });
    this.setState({ ['description']: "" });
    this.setState({ ['foreground_data']: "" });
    this.setState({ ['foregroundformat']: "" });
    this.setState({ ['foreground_filename']: "" });
    this.setState({ ['contextformat']: "2" });
    this.setState({ ['context_data']: "" });
    this.setState({ ['context_filename']: "" });
    this.setState({ ['p_value_cutoff']: 0.001 });
    this.setState({ ['contextformat']: "" });
    this.setState({ ['position_specific']: true });
    this.setState({ ['minimum_occurrences']: 20 });
    this.setState({ ['fold_change_cutoff']: 1.0 });
    this.setState({ ['max_depth']: "" });
    this.setState({ ['width']: 15 });
    this.setState({ ['center_sequences']: false });
    this.setState({ ['multiple_testing_correction']: true });
    this.setState({ ['positional_weighting']: true });
    this.setState({ ['compound_residues']: true });
    this.setState({ ['compound_residue_decomposition']: true });
    this.setState({ ['require_context_id']: true });
    this.setState({ ['extension_direction']: "1" });
    this.setState({ ['redundancylevel']: "1" });
    this.setState({ ['first_protein_only']: true });
    this.setState({ ['originalrowmerge']: "1" });

    document.getElementById('advanced-options-header').innerHTML = 'Advanced options &#9650';
    advancedOptions.style.display = "none";

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
      require_context_id,
      extension_direction,
      redundancylevel,
      first_protein_only,
      originalrowmerge
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
                  checked={foregroundformat == 1}
                  required
                />
                Prealigned text file (<a href='/patterndetection/textfile' download>Example</a>)
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
                  checked={foregroundformat == 3}
                />
                Text file peptide list (<a href='/patterndetection/peptidelist' download>Example</a>)
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
                  checked={contextformat == 2}
                  required
                />
                Swiss-Prot Human
              </label>
            </div>
            <div>
              <label>
                <input
                  type="radio"
                  name="contextformat"
                  id="swissprot-mouse"
                  onChange={this.handleChange}
                  value="3"
                  checked={contextformat == 3}
                />
                Swiss-Prot Mouse
              </label>
            </div>
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
                  + " If pre-aligned sequences are supplied, each supplied sequence MUST be of this length.\n"} />
              <div className="control">
                <input
                  className="input"
                  type="number"
                  name="width"
                  onChange={this.handleChange}
                  value={width}
                  min="1"
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
                    "For example:\n\nCentered:"
                    + "p4-p3-p2-p1-p1'-p2'-p3'-p4\nNon-centered: p1-p2-p3-p4-p5-p6-p7-p8-p9-p10-p11-p12-p13-p14-p15\n"} />
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
                    + " an analysis. When disabled, background frequency are averaged across all posiitons of the"
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
                + " in order to be considered significantly enriched.\n"} />
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
            <div className="field">
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
            </div>
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
                    "Peptides may match more than one position in the context proteome. This setting species how to handle those multiple matches.\n\n"
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