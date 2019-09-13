import React, { Component } from "react";
import PropTypes from "prop-types";
import key from "weak-key";
import getCookie from "./utils";


function formatDateTime(datetime) {
	return (datetime.slice(0, 9) + " " + datetime.slice(11, 19));
}


class JobStatus extends Component {
	static propTypes = {
    	endpoint: PropTypes.string.isRequired,
  	};

	constructor(props) {
		super(props);
		this.state = {
			data: [],
			loaded: false,
			placeholder: 'Loading...',
		};
		this.getJobs = this.getJobs.bind(this);
		this.csrftoken = getCookie('csrftoken');
	}

	componentDidMount() {
		this.getJobs();
	}

	getJobs() {
		fetch(this.props.endpoint, {headers:
			new Headers({ "X-CSRFToken": this.csrftoken})})
    		.then(response => {
        		if (response.status !== 200) {
          			return this.setState({ placeholder: "Something went wrong" });
        		}
        		return response.json();
      		})
      		.then(data => this.setState({ data: data, loaded: true }));
		setTimeout(this.getJobs, 30000);
	}

	render() {
		return  (
			<section>
				<h2 className="subtitle">Recent Jobs</h2>
				<table className="table is-striped">
					<thead>
						<tr>
							<th>Job</th>
							<th>Submitted</th>
							<th>Completed</th>
						</tr>
					</thead>
					<tbody>
						{this.state.data.map((el, index) => (
							<tr key={el.jobcode}>
								<td>{index + 1}</td>
								<td>{formatDateTime(el.submitted)}</td>
								<td>{el.completed ? formatDateTime(el.completed) : "Job queued."}</td>
							</tr>
						))}
					</tbody>
				</table>
			</section>
		)
	}
}

export default JobStatus;