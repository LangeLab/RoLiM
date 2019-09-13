import React, { Component } from "react";
import PropTypes from "prop-types";

class Form extends Component {
  static propTypes = {
    endpoint: PropTypes.string.isRequired
  };
  
  constructor (props) {

  }
  
  state = {
    background_format: ""
  };
  
  handleChange = e => {
    this.setState({ [e.target.name]: e.target.value });
  };
  
  handleSubmit = e => {
    e.preventDefault();
    const { background_format } = this.state;
    const backgroundformat = { background_format };
    const conf = {
      method: "post",
      body: JSON.stringify(backgroundformat),
      headers: new Headers({ "Content-Type": "application/json" })
    };
    fetch(this.props.endpoint, conf).then(response => console.log(response));
  };
  
  render() {
    const { background_format } = this.state;
    return (
      <div className="column">
        <form onSubmit={this.handleSubmit}>
          <div className="field">
            <label className="label">Background Format</label>
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

export default Form;