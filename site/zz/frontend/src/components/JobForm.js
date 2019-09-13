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
    dataupload: "",
    backgroundupload: "",
    alpha: ""
  };

  handleChange = e => {
    this.setState({ [e.target.name]: e.target.value });
  };

  handleUpload = e => {
    this.setState({ [e.target.name]: e.target.files[0] });
    document.getElementById(e.target.name + "-text").textContent = "File selected.";
  };
  
  handleSubmit = e => {
    e.preventDefault();
    let formData = new FormData();
    formData.append('dataupload', this.state.dataupload);
    formData.append('backgroundupload', this.state.backgroundupload);
    formData.append('title', this.state.title);
    formData.append('description', this.state.description);
    formData.append('email', this.state.email);
    formData.append('alpha', this.state.alpha);
    const { title, description, email, dataupload, backgroundupload, alpha } = this.state;
    const job = { title, description, email, dataupload, backgroundupload, alpha };
    const csrftoken = getCookie('csrftoken');
    const conf = {
      method: "post",
      body: formData,
      headers: new Headers({ "X-CSRFToken": csrftoken })
    };
    fetch(this.props.endpoint, conf).then(response => console.log(response));
  };
  
  render() {
    const { title, description, email, dataupload, backgroundupload, alpha } = this.state;
    return (
      <div className="column">
        <form onSubmit={this.handleSubmit} encType="multipart/form-data">
          <div className="field">
            <label className="label">Title</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="title"
                onChange={this.handleChange}
                value={title}
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
            <label className="label">Email</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="email"
                onChange={this.handleChange}
                value={email}
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Select input data file to upload (pre-aligned fixed length sequences, FASTA peptides, or MaxQuant evidence.txt).</label>
            <label className="file-label">
              <div className="control">
                <input
                  className="file-input"
                  type="file"
                  name="dataupload"
                  onChange={this.handleUpload}
                  autoComplete="off"
                />
                <span className="file-cta">
                  <span className="file-icon">
                    <i className="fas fa-upload"></i>
                  </span>
                  <span className="file-label" id="dataupload-text">
                    Choose a file…
                  </span>
                </span>
              </div>
            </label>
          </div>
          <div className="field">
            <label className="label">Select background FASTA file to upload.</label>
            <label className="file-label">
              <div className="control">
                <input
                  className="file-input"
                  type="file"
                  name="backgroundupload"
                  onChange={this.handleUpload}
                  autoComplete="off"
                />
                <span className="file-cta">
                  <span className="file-icon">
                    <i className="fas fa-upload"></i>
                  </span>
                  <span className="file-label" id="backgroundupload-text">
                    Choose a file…
                  </span>
                </span>
              </div>
            </label>
          </div>
          <div className="field">
            <label className="label">Enter p-value threshold for significance test (optional, default=0.01 corrected)</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="alpha"
                onChange={this.handleChange}
                value={alpha}
              />
            </div>
          </div>
          {/*
          <div className="field">
            <label className="label">Select input data file format</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="dataformat"
                onChange={this.handleChange}
                value={dataformat}
                required
              />
            </div>
          </div>
          
          <div className="field">
            <label className="label">Expand peptides in non-prime direction?</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="expand_peptides"
                onChange={this.handleChange}
                value={expand_peptides}
                required
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Include compound residues?</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="compound_residues"
                onChange={this.handleChange}
                value={compound_residues}
                required
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Select reference proteome FASTA file (optional)</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="proteomeupload"
                onChange={this.handleChange}
                value={proteomeupload}
                required
              />
            </div>
          </div>
          <div className="field">
            <label className="label">Select amino acid background frequency file format.</label>
            <div className="control">
              <input
                className="input"
                type="text"
                name="background_format"
                onChange={this.handleChange}
                value={background_format}
                required
              />
            </div>
          </div> */}
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