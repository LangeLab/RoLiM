import React, { Component } from "react";
import PropTypes from "prop-types";
import getCookie from "./utils";

class JobForm extends Component {
  static propTypes = {
    endpoint: PropTypes.string.isRequired
  };
  
  state = {
    title: "",
    description: "",
    email: "",
    foreground_data: "",
    foregroundformat:"1",
    context_data: "",
    p_value_cutoff: "",
    contextformat: "",
    minimum_occurrences: "",
    fold_change_cutoff: "",
    max_depth: "",
    extend_sequences: false,
    width: "",
    center_sequences: true,
    multiple_testing_correction: true,
    positional_weighting: true,
    compound_residues: true,
    compound_residue_decomposition: true
  };

  handleChange = e => {
    this.setState({ [e.target.name]: e.target.value });
  };

  handleCheckboxChange = e => {
    this.setState({ [e.target.name]: e.target.checked });
  }

  handleUpload = e => {
    this.setState({ [e.target.name]: e.target.files[0] });
    document.getElementById(e.target.name + "-text").textContent = "File selected.";
  };
  
  handleSubmit = e => {
    e.preventDefault();
    let formData = new FormData();
    formData.append('foreground_data', this.state.foreground_data);
    formData.append('foregroundformat', this.state.foregroundformat);
    formData.append('context_data', this.state.context_data);
    formData.append('title', this.state.title);
    formData.append('description', this.state.description);
    formData.append('email', this.state.email);
    formData.append('p_value_cutoff', this.state.p_value_cutoff);
    formData.append('minimum_occurrences', this.state.minimum_occurrences);
    formData.append('fold_change_cutoff', this.state.fold_change_cutoff);
    formData.append('extend_sequences', this.state.extend_sequences);
    formData.append('width', this.state.width);
    formData.append('center_sequences', this.state.center_sequences);
    formData.append('multiple_testing_correction', this.state.multiple_testing_correction);
    formData.append('positional_weighting', this.state.positional_weighting);
    formData.append('compound_residues', this.state.compound_residues);
    formData.append('compound_residue_decomposition', this.state.compound_residue_decomposition);

    const {
      title,
      description,
      email,
      foreground_data,
      foregroundformat,
      context_data,
      p_value_cutoff,
      minimum_occurrences,
      fold_change_cutoff,
      extend_sequences,
      width,
      center_sequences,
      multiple_testing_correction,
      positional_weighting,
      compound_residues,
      compound_residue_decomposition
    } = this.state;
    const job = {
      title,
      description,
      email,
      foreground_data,
      foregroundformat,
      context_data,
      p_value_cutoff,
      minimum_occurrences,
      fold_change_cutoff,
      extend_sequences,
      width,
      center_sequences,
      multiple_testing_correction,
      positional_weighting,
      compound_residues,
      compound_residue_decomposition
    };
    const csrftoken = getCookie('csrftoken');
    const conf = {
      method: "post",
      body: formData,
      headers: new Headers({ "X-CSRFToken": csrftoken })
    };

    fetch(this.props.endpoint, conf).then(response => console.log(response));
  };
  
  render() {
    const {
      title,
      description,
      email,
      foreground_data,
      foregroundformat,
      context_data,
      p_value_cutoff,
      minimum_occurrences,
      fold_change_cutoff,
      extend_sequences,
      width,
      center_sequences,
      multiple_testing_correction,
      positional_weighting,
      compound_residues,
      compound_residue_decomposition
    } = this.state;
    return (
      <div className="column">
        <h1>Submit a new job for analysis.</h1>
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
            <label className="label">Title</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="title"
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
            <label className="label">Select input data file to upload.</label>
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
                  onChange={this.handleChange}
                  value="1"
                />
                Prealigned text file
              </label>
            </div>
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
            <div>
              <label>
                <input
                  type="radio"
                  name="foregroundformat"
                  onChange={this.handleChange}
                  value="3"
                />
                Text file peptide list
              </label>
            </div>
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
          </div>

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
          <h2>Advanced options</h2>
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
            <label className="label">Center sequences?</label>
            <div className="control">
              <input
                type="checkbox"
                name="center_sequences"
                onChange={this.handleCheckboxChange}
                value={center_sequences}
                checked={center_sequences}
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Expand and align peptides using context data set?</label>
            <div className="control">
              <input
                type="checkbox"
                name="extend_sequences"
                onChange={this.handleCheckboxChange}
                value={extend_sequences}
                checked={extend_sequences}
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
                onChange={this.handleCheckboxChange}
                value={width}
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Detect compound residue groups?</label>
            <div className="control">
              <input
                type="checkbox"
                name="compound_residues"
                onChange={this.handleCheckboxChange}
                value={compound_residues}
                checked={compound_residues}
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Enable compound residue decomposition?</label>
            <div className="control">
              <input
                type="checkbox"
                name="proteomeupload"
                onChange={this.handleCheckboxChange}
                value={compound_residue_decomposition}
                checked={compound_residue_decomposition}
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Enable multiple testing correction?</label>
            <div className="control">
              <input
                type="checkbox"
                name="multiple_testing_correction"
                onChange={this.handleCheckboxChange}
                value={multiple_testing_correction}
                checked={multiple_testing_correction}
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Enable positional weighting?</label>
            <div className="control">
              <input
                type="checkbox"
                name="positional_weighting"
                onChange={this.handleCheckboxChange}
                value={positional_weighting}
                checked={positional_weighting}
              />
            </div>
          </div>
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