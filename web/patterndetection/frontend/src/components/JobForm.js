import React, { Component } from "react";
import PropTypes from "prop-types";
import getCookie from "./utils";

function HelpText(props) {
  render() {
    return (
      <button onClick={() => alert(props.text)}>?</button>
    );
  }
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
    foregroundformat: "1",
    foreground_filename: "",
    context_data: "",
    contextformat: "",
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
  };

  handleChange = e => {
    this.setState({ [e.target.name]: e.target.value });
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

    const {
      title,
      description,
      email,
      foreground_data,
      foregroundformat,
      foreground_filename,
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
      require_context_id
    } = this.state;
    const job = {
      title,
      description,
      email,
      foreground_data,
      foregroundformat,
      foreground_filename,
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
      require_context_id
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
    this.setState({ ['foregroundformat']: 1 });
    this.setState({ ['foreground_filename']: "" });
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
      require_context_id
    } = this.state;
    return (
      <div className="column">
        <h4 className="title is-4">Submit a new job for analysis.</h4>
        <form onSubmit={this.handleSubmit} encType="multipart/form-data">
          <div className="field">
            <label className="label">Email</label>
            <div className="control">
              <input
                className="input"
                type="text"
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
                  required
                />
                Prealigned text file (<a href='/patterndetection/textfile' download>Example</a>)
              </label>
            </div>
            {/*
            <div>
              <label>
                <input
                  type="radio"
                  name="foregroundformat"
                  onChange={this.handleChange}
                  value="2"
                />
                Prealigned FASTA file
              </label>
            </div>
            */}
            <div>
              <label>
                <input
                  type="radio"
                  name="foregroundformat"
                  id="peptidelist"
                  onChange={this.handleChange}
                  value="3"
                />
                Text file peptide list (<a href='/patterndetection/peptidelist' download>Example</a>)
              </label>
            </div>
            {/*
            <div>
              <label>
                <input
                  type="radio"
                  name="foregroundformat"
                  onChange={this.handleChange}
                  value="4"
                />
                FASTA peptide list
              </label>
            </div>
            <div>
              <label>
                <input
                  type="radio"
                  name="foregroundformat"
                  onChange={this.handleChange}
                  value="6"
                />
                MaxQuant "evidence.txt" file
              </label>
            </div>
            */}
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
                  Center sequence position numbers? <HelpText text="For example:\n\nCentered: p4-p3-p2-p1-p1'-p2'-p3'-p4\nNon-centered: p1-p2-p3-p4-p5-p6-p7-p8" />
                </label>
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
              </div>
            </div>
            <br />
            <div className="field">
              <label className="label">P-value threshold.</label>
              <div className="control">
                <input
                  className="input"
                  type="text"
                  name="p_value_cutoff"
                  onChange={this.handleChange}
                  value={p_value_cutoff}
                />
              </div>
            </div>
            <div className="field">
              <label className="label">Minimum occurrences</label>
              <div className="control">
                <input
                  className="input"
                  type="text"
                  name="minimum_occurrences"
                  onChange={this.handleChange}
                  value={minimum_occurrences}
                />
              </div>
            </div>
            <div className="field">
              <label className="label">Fold change cutoff</label>
              <div className="control">
                <input
                  className="input"
                  type="text"
                  name="fold_change_cutoff"
                  onChange={this.handleChange}
                  value={fold_change_cutoff}
                />
              </div>
            </div>
            <div className="field">
              <label className="label">Enter desired width of expanded sequences.</label>
              <div className="control">
                <input
                  className="input"
                  type="text"
                  name="width"
                  onChange={this.handleChange}
                  value={width}
                />
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
      </div>
    );
  }
}

export default JobForm;